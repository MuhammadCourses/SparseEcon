
## This files saves details on how to get started with using Sparse Econ repository by Schaab and Zhang. In particular, how to setp up adaptative vs uniform grids. 

# Grid Generator (`setup_grid.m`)

`setup_grid.m` builds either a **dense** (full tensor) or **sparse** (Smolyak) grid for a \(d\)-dimensional state space.

---

## Quick Examples (Huggett model)

```matlab
% Dense 65×65 grid over (a,z):
G_dense = setup_grid(0, [6 6], param.min, param.max, ...
                     'NamedDims',{1 2}, 'Names',{'a','z'}, 'DxxDims',2);
G_dense.dx = G_dense.da * G_dense.dz;        % area element

% Sparse Smolyak grid (level 6):
G = setup_grid(6, [0 0], param.min, param.max, ...
               'NamedDims',{1 2}, 'Names',{'a','z'}, 'DxxDims',2);
```

---
## How It Works

1. **1-D nested meshes**  
   Level \($\ell$\) has \($2^{\ell}+1$\) equally spaced points on \([0,1]\).

2. **Dense grid** (`n=0`, large `surplus`)  
   Keeps **all** level combinations up to \($\ell_i\le L$\) → \(($2^L+1)^d$\) nodes.

3. **Sparse grid (Smolyak)**  
   Keeps multi-indices \($\ell$\) satisfying  
   $$
     \sum_i \max\{\ell_i-\text{surplus}(i),0\}\le n,
   $$ 
   giving polynomial (not exponential) growth in node count.

4. **Scaling**  
   Each coordinate is mapped to \([$\text{min}_i,\text{max}_i$]\).

5. **Extras**  
   Optional tags (`'NamedDims','Names'`) and derivative matrices (`'DxxDims','DxyDims'`).

---

## Input Summary

| Argument     | Meaning                                                                                         |
|--------------|-------------------------------------------------------------------------------------------------|
| **`n`**      | Smolyak level (higher = finer sparse grid).                                                     |
| **`surplus`**| 1×d vector of per-dimension boosts; `surplus(i)>0` forces full tensor resolution in dim i.      |
| **`min`**    | 1×d vector of lower bounds in economic units.                                                   |
| **`max`**    | 1×d vector of upper bounds.                                                                     |
| **Options**  | `'NamedDims','Names'` – label dims; `'DxxDims','DxyDims'` – pre-compute 2nd / mixed derivatives.|

---

## Typical Choices

| Goal                          | Example settings                              |
|-------------------------------|-----------------------------------------------|
| Sparse, 2–3 dims              | `n = 6;  surplus = zeros(1,d);`               |
| Extra accuracy in assets dim  | `n = 6;  surplus = [2 0 … 0];`                |
| Dense 65×65 grid              | `n = 0;  surplus = [6 6];`                    |

---

## Projection Matrix: Sparse → Dense Projection
```matlab
G.BH_dense = get_projection_matrix(G_dense.grid, G_dense.lvl, G);
% where G_dense.grid is grid of dense points, G_dense_lvl is level of sparse grid that was used in each dimension. It comes automatically from `gen_sparse_grid.`
```
BH_dense is a sparse matrix that maps any function stored on the sparse grid G to the dense grid G_dense:
```matlab
f_dense = G.BH_dense * f_sparse;   % one mat-vec, no loops
```

Once we know how to create grids, their main use is going to be in taking derivatives. Here are how to do that:
## Derivative Operators


## 4. Finite‑Difference Operator to take derivative of functions of diffusions (`FD_operator.m`)

`FD_operator` assembles the semi‑discrete drift–diffusion operator

$$
L f \;=\;
\sum_{k\in\texttt{dims}} \mu_k\,\partial_{x_k} f
\;+\;
\frac12 \sum_{k\in\texttt{dims}} \sigma_k^{2}\,\partial_{x_k x_k}^{\,2} f
\;+\;
\sum_{k<\ell} \sigma_k\sigma_\ell\,\partial_{x_k x_\ell}^{\,2} f ,
$$

on the grid `G`.

---

### Call

```matlab
[ A , b ] = FD_operator(G, mu, sigma, dims, BC_name)
```

| Argument   | Size / Type | Description |
|------------|-------------|-------------|
| `G`        | struct      | Grid with pre‑built FD stencils (`D1F`,`D1B`,`D2`,`D22`). |
| `mu`       | `J × d'`    | Point‑by‑point **drift** along axes `dims`. |
| `sigma`    | `J × d'`    | Point‑by‑point **volatilities** (standard deviation, not variance). |
| `dims`     | vector      | Subset of grid dimensions the columns of `mu`,`sigma` refer to. |
| `BC_name`  | string *(opt.)* | Which boundary‑condition set in `G` to use (`'main'` default). |

*Returns*

| Output | Meaning |
|--------|---------|
| `A`    | `J × J` sparse (or dense) matrix for interior nodes. |
| `b`    | `J × 1` constant vector that absorbs boundary terms. |

---

### Assembly Logic

| Part | Action |
|------|--------|
| **First derivatives** | Upwind automatically: forward stencil if `mu ≥ 0`, backward if `mu < 0`. |
| **Second derivatives** | Adds $$\tfrac12\sigma^2 D2$$. Requires that axis be listed in `DxxDims`. |
| **Mixed derivatives**  | Adds $$\sigma_k\sigma_\ell D22$$. Requires `[k ℓ]` in `DxyDims`. |
| **Sparse vs dense**    | Chooses sparse or dense storage to match `G.sparse`. |

---

### Derivatives when process has discrete types or is poisson


### Role of the Constant Vector `b`

Finite‑difference stencils next to the boundary reference **ghost points**.  
Because ghost values are **fixed by the boundary condition**, their contribution is moved to the right‑hand side:

$$
A\,f_{\text{unknown}} \;+\; b \;=\; \text{RHS}.
$$

| BC type  | Ghost‑point replacement | Contribution to `b` |
|----------|-------------------------|---------------------|
| Dirichlet $f = f_{\text{BC}}$ | Substitute fixed value in stencil | Adds $\pm c\,f_{\text{BC}}$ |
| Neumann $\partial_x f = g_{\text{BC}}$ | Ghost value $= f_1 \pm h\,g_{\text{BC}}$ | Adds $\pm c\,h\,g_{\text{BC}}$ |

Hence you **solve**

```matlab
A * f_interior = rhs - b;
```

without rebuilding `A` when boundary data change.

---

### Shortcuts

* **1‑D, drift‑only**

  ```matlab
  A = FD_operator(G , mu);   % upwind drift matrix only
  ```

* **Stationary HJB / Fokker–Planck**

  ```matlab
  f = A \ (rhs - b);        % direct solve
  ```

---

### Workflow Summary

1. Generate **dense** or **sparse** grid with `setup_grid`.  
2. Evaluate **drift** (`mu`) and **volatility** (`sigma`) at grid nodes.  
3. Call `FD_operator` to obtain matrices `A, b`.  
4. Solve or time‑step your PDE:  
   * steady state → `f = A \ (rhs - b)`  
   * dynamics     → implicit / explicit schemes with `A`.

See inline comments in `FD_operator.m` for implementation details.

## Taking Derivative of simple analytical function `deriv_sparse.m`

Compute finite-difference derivatives on a Smolyak sparse grid.

## Usage

```matlab
deriv = deriv_sparse(G, f, k, operator, name)
```
| Argument | Description |
|----------|-------------|
| `G` | Grid struct from setup_grid (sparse grid). |
| `f` | J×1 vector of values at grid nodes. |
| `k` | Scalar index for univariate operators (1, 2, …).<br>Two-element vector [i j] for mixed derivative D22. |
| `operator` | `'D1F'` forward first derivative<br>`'D1B'` backward first derivative<br>`'D2'` second derivative<br>`'D22'` mixed second derivative |
| `name` | (optional): Boundary-condition set (default 'main'). |

**Returns**

| Output | Description |
|--------|-------------|
| `deriv` | J×1 vector of derivative values at each grid node. |

**Example Usage**
```matlab
VaF = deriv_sparse(G, V, 1, 'D1F');
VaB = deriv_sparse(G, V, 1, 'D1B');
```