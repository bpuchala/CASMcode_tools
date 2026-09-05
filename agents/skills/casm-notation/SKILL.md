---
name: casm-notation
description: CASM symmetry/group-theory notation, terminology, and math-symbol conventions (little group, prim/supercell factor group, k-points, factor systems, etc.) for reading, writing, and discussing CASM docs and code. Use whenever writing docstrings/comments or explaining symmetry, irrep, or group-theory results in a CASM package.
---

# CASM notation and terminology

Vocabulary, notation, and math-symbol conventions for CASM symmetry and
group-theory work. For *how* to write docs (voice, plain language, phrasing),
see the `casm-dev` skill's `documentation-style.md`.

Distinguish the exact (often infinite) mathematical group from the finite object
CASM stores to represent it. CASM typically stores a quotient/factor group as
its finite set of coset representatives; make clear which one is meant.

## Preferred vocabulary

Use conventions from crystallography (Bradley & Cracknell) rather than pure
mathematics:

- 2-cocycle -> factor system
- the group of cocycles mod coboundaries -> cohomology group
- group of phases / multipliers -> multiplier group
- projective representation -> projective (ray) representation, or ray representation

Group theory terms that we expect will be understood (or will provide background for):

- "group", "subgroup"
- "coset", "coset representative"
- "multiplication table"
- "translation group"
- "space group", "point group"
- "quotient group" / "factor group" (synonyms)
- "representation", "matrix representation"
- "invariant subspace"
- "reciprocal lattice", "wave vector", "k-point"
- "irreducible representation" / "irrep"
- "character" (of a representation)
- "symmetry-adapted basis"
- "orbit" (of symmetry-equivalent k-points, modes, sites, clusters, etc.)
- "wedge"
- "Reynolds operator"
- "equivalence map"

Terms that may cause AI tools confusion:

- "invariant subgroup" -- in the CASM API (e.g. `make_invariant_subgroup`) and documentation, this is the stabilizer of an orbit element: the subgroup of operations that leave a particular object (site, cluster, k-point, configuration) invariant. NOTE: this is NOT the classical / Bradley & Cracknell meaning, where "invariant subgroup" = normal (self-conjugate) subgroup.

Mathematical terms to avoid and preferred form:

- "cohomologous" / "differ by a coboundary" / "gauge-equivalent" -> "equivalent up to a choice of basis or phase factor", or just "give the same irreps"
- "gauge" / "gauge transformation" -> "a choice of phase factor"
- "chain", "cochain", "cycle", "cocycle" -> avoid; describe the factor system or multiplier group directly instead
- "TRIM" / "time-reversal invariant momentum" -> "a k-point where k = -k modulo a reciprocal lattice vector"; for the negation, "a k-point where k != -k" rather than "non-TRIM"
- "stabilizer" -> "invariant subgroup"

Terms that may be used, but must be explained:

- "isotypic component"
- "factor system"
- "multiplier group"
- "projective representation"
- "corepresentation" -- a representation of a group that includes time-reversal
  (antiunitary) operations

Preferred forms (spelling / usage):

- "k-point" in prose ("kpoint" only as a code identifier)
- "symmetry-adapted" (hyphenated), not "symmetry adapted"
- "time-reversal" (adjective) / "time reversal" (noun); "time_reversal" only as a
  code identifier
- introduce "irreducible representation (irrep)" on first use, then "irrep"
- "factor group" instead of "quotient group" when referring to the "prim factor group" or "supercell factor group"

## Notation and CASM terminology

- $G$: The infinite space group of the crystal and allowed variable space.
- $T$: The infinite translation group of the crystal.
- $G/T$: The finite quotient group of the crystal. In CASM, this is the *prim factor group*. 
- $G_S$: The infinite space group of the user-specified supercell, $S$. This is the subgroup of $G$ containing all operations that leave the supercell lattice vectors invariant. It does include operations with fractional translations of the primitive cell (glide/screw).
- $G_S/T$: The finite quotient group of supercell, $S$, with respect to the translation group of the crystal. This is the subgroup of $G/T$ containing only operations that leave the supercell lattice invariant. In CASM, this is the *supercell factor group*.
- $T_S$: The infinite translation group of the supercell $S$.
- $G_S/T_S$: The finite quotient group of supercell, $S$, with respect to the translation group of the supercell. This subdivides the cosets of $G_S/T$ into $N^{u}$ smaller cosets, one per translation inside the supercell. In CASM, these are represented by *SupercellSymOp*.
- $G^{\mathbf{k}}$: The infinite little group, the subgroup of $G$ which leaves $\mathbf{k}$ invariant ($\mathbf{k} = R\mathbf{k}$ modulo the prim reciprocal lattice vectors). It does include operations with fractional translations of the primitive cell (glide/screw). In CASM, this is the *little group*.
- $T^{\mathbf{k}}$: The infinite translation group of $\mathbf{k}$, a subgroup of $T$, it consists of translations $\mathbf{t}$ for which the Bloch phase factor is unity ($T^{\mathbf{k}} = \{ \{E|\mathbf{t}\} \in T \mid e^{-i\mathbf{k}\cdot\mathbf{t}} = 1 \}$). $T^{\mathbf{k}}$ is the infinite translation group of the smallest supercell commensurate with $\mathbf{k}$, which we will call $S^{\mathbf{k}}$.
- $G^{\mathbf{k}}/T^{\mathbf{k}}$: The finite quotient group of the little group of $\mathbf{k}$ and the translation group of $\mathbf{k}$.
- $\bar{G}^{\mathbf{k}}$: The *little co-group*, the finite point group of $\mathbf{k}$, consisting of the rotational parts $R$ of the little group operations $\{R|\tau\} \in G^{\mathbf{k}}$ (those with $R\mathbf{k} = \mathbf{k}$ modulo the prim reciprocal lattice vectors). Unlike the little group $G^{\mathbf{k}}$, it discards the translations. Keep distinct from the little group.
- $G_S^{\mathbf{k}}$: The infinite supercell little group, this is the group formed by the intersection of the little group and the space group of the user-specified supercell, $G_S^{\mathbf{k}} = G^{\mathbf{k}} \cap G_S$. These are space group operations that leave both $\mathbf{k}$ and the supercell lattice vectors invariant.
- $G_S^{\mathbf{k}} / T_S$: The finite quotient group of the supercell little group and the supercell translation group. In CASM, this is the *supercell little group*.

## Math symbols

Use these consistently in docstrings, comments, and LaTeX/`:math:`. These match existing usage; only adopt new forms going forward (do not churn existing docs solely to convert notation).

- Vector: bold lowercase, e.g. $\mathbf{v}$, $\mathbf{k}$, $\mathbf{x}$, $\mathbf{t}$.
- Matrix: plain (italic) uppercase, e.g. $M$, $D$, $Q$, $R$.
- Unit vector: a hat over the vector symbol, e.g. $\hat{\mathbf{e}}$, $\hat{\mathbf{n}}$.
- Complex conjugate: superscript star, e.g. $M^*$, $v^*$. Reserve $^*$ for complex conjugation only -- use a prime ($t'$) or hat ($\hat{L}_k$) for chosen/distinguished objects or group extensions, never $^*$, and do not use an overbar (the overbar is the little co-group, e.g. $\bar{G}^{\mathbf{k}}$).
- Transpose: superscript $\top$, e.g. $M^\top$ (LaTeX `\top`). The superscript is distinct from the transformation matrix $T$ in $S = L T$. (Older docs use `^T`; leave them, use `\top` going forward.)
- Conjugate transpose (adjoint): a dagger, e.g. $M^\dagger$ in LaTeX, `†` in plain-text comments (e.g. `Q^†`). Never $^*$ or $^H$. In code: Eigen `.adjoint()`, NumPy `.conj().T`.
- Inverse: $M^{-1}$.
- Pseudoinverse (Moore-Penrose): $M^{+}$. In code: Eigen `.completeOrthogonalDecomposition().pseudoInverse()`, NumPy `np.linalg.pinv`. 
