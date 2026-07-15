# Multipole Interactions: Automated Derivation of Cartesian Multipole Terms

[![License: CC BY-NC 4.0](https://img.shields.io/badge/License-CC%20BY--NC%204.0-blue.svg)](https://creativecommons.org/licenses/by-nc/4.0)

> **Author**: [Judith M. Leson](https://orcid.org/0009-0002-5627-2673)  
> **Repository**: [https://github.com/jmleson/multipole_interactions](https://github.com/jmleson/multipole_interactions)  
> **Related Work**: Dissertation: *A Quantum-Chemical Analysis of Long-Range Dimer Interactions Arising From Triplet Excited States of Monocyclic Aromatics* (Judith M. Leson, University of Duisburg-Essen, 2026)  
> **Data DOI**: [10.71955/DUEDATA-2026-MR4YY63J](https://doi.org/10.71955/DUEDATA-2026-MR4YY63J)  
> **Reference**: [Stone, A. (2013). *The Theory of Intermolecular Forces*, 2nd ed., Oxford University Press.](https://doi.org/10.1093/acprof:oso/9780199672394.001.0001)

---

## 📌 Overview

This repository provides automated Python code to symbolically derive Cartesian multipole interaction terms for two molecules in a stacked geometry along the Z-axis. 
The focus is on generating **traceable, step-by-step derivations** of interaction formulas up to high multipole orders.

The code is intended for **computational and theoretical chemists** who need to work with multipole expansions of different orders and want to trace, verify or adapt the underlying equations.

[//]: # (The code can be of use for researchers in theoretical and computational chemistry who need traceable and adaptable expressions for multipole interactions.)

---
## ✅ Key Features

- ✅ **Symbolic derivation** of Cartesian multipole tensors from $ \nabla_\alpha \nabla_\beta \cdots (1/R) $
- ✅ **Stepwise simplification** for stacked systems (intermolecular distance along Z)
- ✅ **Support for arbitrary multipole orders** (dipole, quadrupole, octopole, hexadecapole, etc.)
- ✅ **Automated LaTeX File generation** for tensor definitions and interaction terms
- ✅ **Human-readable derivations in compiled PDFs** showing the intermediate steps

---
## 🛠️ How It Works
### Theoretical Foundation
The derivation is based on the general definition of Cartesian multipole tensors:  
$$T_{\alpha\beta\cdots\nu} = \nabla_\alpha \nabla_\beta \cdots \nabla_\nu \left( \frac{1}{R} \right) \quad, $$  
with $$R = |\vec{R}|=\sqrt{R_x^2 + R_y^2 + R_z^2} \quad .$$  

These tensors $T$ can be used to compute the interaction energy $E_{interaction}$ between two molecules A and B via a series expansion (Stone 2013). 
Considering systems without permanent charges, the multipole expansion takes the form:   
$$E_{interaction} = $$
$$T_{\alpha\beta} \cdot (-\mu_\alpha^B \mu_\beta^A)$$   
$$+ T_{\alpha\beta\gamma} \cdot \left( +\mu_\alpha^B \frac{1}{3} \Theta_{\beta\gamma}^A - \frac{1}{3} \Theta_{\alpha\beta}^B \mu_\gamma^A \right) $$  
$$+ T_{\alpha\beta\gamma\delta} \cdot \left( -\frac{1}{15} \mu_\alpha^B \Omega_{\beta\gamma\delta}^A + \frac{1}{9} \Theta_{\alpha\beta}^B \Theta_{\gamma\delta}^A - \frac{1}{15} \Omega_{\alpha\beta\gamma}^B \mu_\delta^A \right) $$  
$$+ T_{\alpha\beta\gamma\delta\epsilon} \cdot \left( +\frac{1}{105} \mu_\alpha^B \Phi_{\beta\gamma\delta\epsilon}^A - \frac{1}{45} \Theta_{\alpha\beta}^B \Omega_{\gamma\delta\epsilon}^A + \frac{1}{45} \Omega_{\alpha\beta\gamma}^B \Theta_{\delta\epsilon}^A - \frac{1}{105} \Phi_{\alpha\beta\gamma\delta}^B \mu_\epsilon^A \right) $$  
$$+ \cdots$$  


> **Notation**:
> - $\mathbf{T}_{\alpha\beta\cdots}$: Cartesian multipole tensor  
> - $\mu$: dipole moment  
> - $\Theta$: quadrupole moment  
> - $\Omega$: octopole moment  
> - $\Phi$: hexadecapole moment  
> - Superscripts $A$, $B$: molecule labels

Each term in the expansion corresponds to specific multipoles order.
The code derives these terms systematically using symbolic differentiation and by enforcing the traceless conditions of the respective multipole components.


### Geometric Constraints
In this work, we consider stacked systems (e.g., the "sandwich" benzene dimer), where the intermolecular distance is aligned along the $z$-axis. Therefore, the derivations assume:  
$$R_x = R_y = 0 \Rightarrow R = R_z$$  


Additionally, we assume that certain multipole components vanish due to symmetry in systems such as the benzene, chlorobenzene, and pyrazine dimer:  
$$\mu_\alpha = 0 \quad \text{if } \alpha \neq y$$  
$$\Theta_{\alpha\beta} = 0 \quad \forall \alpha \neq \beta$$  
$$\Omega_{\alpha\beta\gamma} \quad \text{if } \alpha\beta\gamma \notin \left\lbrace xxy, yyy, yzz \right\rbrace$$  
$$\Phi_{\alpha\beta\gamma\delta} \quad \text{if } \alpha\beta\gamma\delta \notin \left\lbrace xxxx, xxyy, xxzz, yyyy, yyzz, zzzz \right\rbrace$$

For systems satisfying these geometric and symmetry constraints, the code performs stepwise derivations and simplifications to yield final interaction equations for specific multipole combinations.
We note that while these constraints apply only to certain symmetric systems, the code remains useful for general use: it generates derivations step-by-step, enabling full traceability and verification of used simplifications.


### Example 

[//]: # (In our case, we considered a stacked benzene &#40;"sandwich"&#41;, where the intermolecular distance goes along the Z coordinate. )
[//]: # (For such systems, the code performs stepwise derivations and simplifications, that result in interaction equations for different multipole interactions. )
For instance, the interaction between two dipoles (rank 1) in a stacked benzene dimer ("sandwich" structure) can be calculated as follows: 

$$- 1 \cdot T_{\alpha{}\beta{}} \mu{}^{B}_{\alpha{}} \mu{}^{A}_{\beta{}} = R^{-5} \cdot (-1) \cdot  \left[ \begin{array}{c} 
 		+ 3 \cdot R_{\alpha{}} \cdot R_{\beta{}} \cdot \mu{}^{B}_{\alpha{}} \cdot \mu{}^{A}_{\beta{}} \\
		 - 1 \cdot R^{2} \cdot \delta{}\left(\alpha{}, \beta{}\right) \cdot \mu{}^{B}_{\alpha{}} \cdot \mu{}^{A}_{\beta{}}
 \end{array}\right] 
$$
Simplify deltas:
$$- 1 \cdot T_{\alpha{}\beta{}} \mu{}^{B}_{\alpha{}} \mu{}^{A}_{\beta{}} = R_{z}^{-5} \cdot -1 \cdot  \left[ \begin{array}{c} 
 		+ 3 \cdot R_{\alpha{}} \cdot R_{\beta{}} \cdot \mu{}^{B}_{\alpha{}} \cdot \mu{}^{A}_{\beta{}} 
 \\ 		 - 1 \cdot R_{z}^{2} \cdot \delta{}\left(\alpha{}, \alpha{}\right) \cdot \mu{}^{B}_{\alpha{}} \cdot \mu{}^{A}_{\alpha{}}
 \end{array}\right] $$
$$ \Rightarrow 
 R_{z}^{-5} \cdot -1 \cdot  \left[ \begin{array}{c} 
 		+ 3 \cdot R_{\alpha{}} \cdot R_{\beta{}} \cdot \mu{}^{B}_{\alpha{}} \cdot \mu{}^{A}_{\beta{}} 
 \\ 		 - 1 \cdot R_{z}^{2} \cdot \mu{}^{B}_{\alpha{}} \cdot \mu{}^{A}_{\alpha{}}
 \end{array}\right] 
$$
Simplify R:
$$- 1 \cdot T_{\alpha{}\beta{}} \mu{}^{B}_{\alpha{}} \mu{}^{A}_{\beta{}} = R_{z}^{-5} \cdot -1 \cdot  \left[ \begin{array}{c} 
 		+ 3 \cdot R_{z} \cdot R_{z} \cdot \mu{}^{B}_{z} \cdot \mu{}^{A}_{z} 
 \\ 		 - 1 \cdot R_{z}^{2} \cdot \mu{}^{B}_{\alpha{}} \cdot \mu{}^{A}_{\alpha{}}
 \end{array}\right] $$
$$ \Rightarrow 
 R_{z}^{-5} \cdot -1 \cdot  \left[ \begin{array}{c} 
 		- 1 \cdot R_{z}^{2} \cdot \mu{}^{B}_{\alpha{}} \cdot \mu{}^{A}_{\alpha{}}
 \end{array}\right] $$
Simplify Dipoles:
$$- 1 \cdot T_{\alpha{}\beta{}} \mu{}^{B}_{\alpha{}} \mu{}^{A}_{\beta{}} = R_{z}^{-5} \cdot -1 \cdot  \left[ \begin{array}{c} 
 		- 1 \cdot R_{z}^{2} \cdot \mu{}^{B}_{y} \cdot \mu{}^{A}_{y}
 \end{array}\right] $$
Combining everything:
$$- 1 \cdot T_{\alpha{}\beta{}} \mu{}^{B}_{\alpha{}} \mu{}^{A}_{\beta{}} =  \left[ \begin{array}{c} 
 		+ 1 \cdot R_{z}^{-3} \cdot \mu{}^{B}_{y} \cdot \mu{}^{A}_{y}
 \end{array}\right] 
$$


---

## 🔧 How to Use
### 📂 Structure
- `run.py`: Main script that generates tensor definitions and interaction examples
- `example.ipynb`: Introduction how to use this library
- `requirements.txt`: Needed python libraries
- `src`: Folder for source files
- `run_latex.sh`: Compiles `.tex` files into PDFs (output in `build/`)
- `build`: Folder for LaTeX compilation files

### 📦 Dependencies

This project requires:
- `sympy` 
- `itertools` 
- `collections`

Install with:
```bash
python3 -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
```


### 💻 Running the Code
Exemplary use of the code is given in `run.py`. For a more interactive example we give the jupyter notebook `example.ipynb`.  

In general, the results of a different order for a tensor can be gained by running: 
```python
save_lines_as_latex_file(filename="tensor_definitions.tex", importable_tex=False, content=[Tensor(order=order).to_tex(importable_tex=True)])
```
or for multiple cases as: 
```python
latex_doc = []
for order in range(3):
    t = Tensor(order=order)
    latex_doc.append( t.to_tex(importable_tex=True) )

save_lines_as_latex_file(filename="tensor_definitions.tex", importable_tex=False, content=latex_doc)
``` 

The interaction between two multipoles in a stacked system can be calculated by giving their rank. 
For instance, the function call for the interaction between two dipoles is: 
```python
s = MultipoleInteraction(multipole_order_1 = 1, multipole_order_2 = 1)
s.simplify_in_latex_steps(importable_tex=False)
```
The parameter ```importable_tex``` determines whether the generated LaTeX file includes a document header (e.g., \documentclass, \begin{document}).
- True: Full LaTeX document (ready to compile).
- False: Only the content (e.g., for inclusion in another LaTeX document).
