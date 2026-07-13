# Documentation style

How to write CASM docstrings, comments, and prose. For the vocabulary, notation,
and math symbols themselves (little group, factor group, `G`/`T`, etc.), use the
`casm-notation` skill.

## Voice and altitude

- Prefer plain language, or pair technical language (for precision) with plain
  language for understanding. The reader is a materials scientist/crystallographer
  following Bradley & Cracknell conventions, not a pure mathematician.
- Use crystallography (Bradley & Cracknell) vocabulary rather than pure
  mathematics — see the `casm-notation` skill for the preferred terms and the
  math terms to avoid (e.g. "gauge" -> "a choice of phase factor").
- Distinguish the exact (often infinite) mathematical group from the finite
  object CASM stores to represent it; make clear which one is meant.
- **NO NEWLINES IN MARKDOWN PARAGRAPHS OR LIST ITEMS** so it is easier for a
  human to edit.

## Phrasing

- **Positive phrasing** — describe what the code *does*, not what it avoids,
  doesn't need, or doesn't do. Drop comparisons to removed or alternative
  approaches (e.g. write "built from prim data via `make_projective_dof_space_rep`
  with `phase_convention='supercell'`", not "built without any supercell-scale
  operator" or "unlike vX it doesn't ...").
- **Don't overstate** — describe what the algorithm does operationally (e.g.
  "each block is checked for irreducible invariance; components that fail are
  re-split"). Hedge an unproven "why" with its observed scope ("observed for some
  nonsymmorphic non-Γ orbits, e.g. diamond with a glide"); do not assert a
  hypothesis drawn from a few cases as established fact.
- **No invented terms** — do not coin ad-hoc labels (e.g. "m_true", "action").
  Use names already in the codebase, or name the concrete object directly (e.g.
  "consistent with the supercell DoF matrix representation `D(g)` produced by
  `make_dof_space_rep`", not "supercell-action-consistent"). Before using a noun
  for a concept, check it appears in the codebase or the `casm-notation` skill.

## Editing existing prose

- **Preserve the user's prose.** When editing docstrings/comments the user
  authored, make the smallest delta that accomplishes the change: delete what
  must go, add only what is strictly necessary, and leave their wording intact
  everywhere else. If a broader rewrite would improve flow, propose it explicitly
  rather than folding it silently into the edit.
