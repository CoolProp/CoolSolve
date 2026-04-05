# CoolSolve Documentation Strategy

## Context

The CoolSolve repository already contains structured Markdown documentation:

- `README.md` — project overview, features, build instructions
- `docs/debugging_models.md` — how to diagnose solver failures
- `docs/deployment_ubuntu_apache.md` — server deployment guide
- `docs/gui.md` — GUI architecture reference
- `docs/solver_roadmap.md` — solver algorithm roadmap
- `docs/symbolic_redecomposition.md` — symbolic reduction internals

The goal is to surface this content as official, navigable documentation reachable from:

1. The GitHub repository page
2. The CoolSolve GUI (a link in the toolbar or help panel)
3. Optionally the dedicated server

**Hard constraints:**

- **No duplication** — one source of truth; content lives in `.md` files only
- **Minimal external dependency** — avoid services that could disappear or require accounts
- **Syntax reference** — a new page documenting the EES/CoolSolve language must be added

---

## Option A — GitHub Pages with MkDocs (static site, CI-generated)

The `.md` files are rendered into a static HTML site by [MkDocs](https://www.mkdocs.org/) (or MkDocs-Material) and deployed to the `gh-pages` branch automatically by a GitHub Actions workflow on every push to `main`.

**How it works:**

1. Add a `mkdocs.yml` at the repo root defining the navigation structure.
2. A GitHub Actions workflow (`deploy-docs.yml`) runs `mkdocs build` and deploys to `gh-pages`.
3. GitHub Pages serves the static site at `https://<org>.github.io/CoolSolve/`.
4. The GUI toolbar links to that URL.

**What you maintain:** only the `.md` files. The HTML is regenerated automatically.

**Dependencies introduced:**
- `mkdocs` Python package (only needed in CI, not in the C++ build)
- GitHub Actions (already available on any public GitHub repo)
- The `gh-pages` branch (managed entirely by the workflow, never edited by hand)

---

## Option B — GitHub Pages with Pandoc (no Python framework)

Instead of MkDocs, a `Makefile` or shell script converts the `.md` files to HTML using [Pandoc](https://pandoc.org/) and pushes the result to `gh-pages`. A hand-written `index.html` provides navigation.

**Dependencies introduced:**
- `pandoc` (single binary, available in most CI images)
- A small shell script (~30 lines) to stitch the pages together

**Trade-off vs A:** more manual effort to maintain navigation; no search, no theme — but zero Python and zero framework lock-in.

---

## Option C — GitHub Pages with no build step (plain Markdown rendering)

GitHub natively renders Markdown files in the repository browser. An `index.md` with links to each doc page, combined with the GitHub Pages "Jekyll" processor (which is on by default), produces a minimal styled site with **zero tooling** — no CI workflow, no build step, nothing to install.

**How it works:**

1. Enable GitHub Pages from the `main` branch, `/ (root)` or `/docs` folder.
2. Each `.md` is already a rendered HTML page at its path under `https://<org>.github.io/CoolSolve/`.
3. Add a `_config.yml` (3 lines: title, description, theme name) to pick a Jekyll theme.

**Dependencies introduced:** none. GitHub's server runs a safe subset of Jekyll automatically.

**Limitations:**
- No custom navigation sidebar
- No search
- No syntax highlighting themes beyond GitHub's default
- Jekyll's Markdown flavour occasionally differs from GFM

---

## Option D — Self-hosted: serve the `.md` files via the CoolSolve binary itself

The CoolSolve GUI binary already embeds and serves static assets (the React SPA). The same mechanism (`embed_assets.cmake`) can embed the `.md` files and serve them as HTML through a `/docs/` route, with the browser rendering them using a lightweight JS Markdown renderer (e.g. [marked.js](https://marked.js.org/), a ~50 kB single file, no npm dependency at runtime).

**How it works:**

1. The `.md` files are embedded at build time alongside the GUI assets.
2. A `/docs/{page}` HTTP route in `server.cpp` returns the raw Markdown content.
3. A `docs.html` shell page (also embedded) fetches the requested `.md` via `fetch()` and renders it client-side with `marked.js`.
4. The GUI toolbar "Help / Docs" button opens `http://localhost:8550/docs/`.

**Dependencies introduced:**
- `marked.js` — one minified JS file (~50 kB), bundled in the repo, no CDN call at runtime
- A small HTTP route handler in `server.cpp` (~20 lines)

**No GitHub Pages needed.** Documentation works offline, on the dedicated server, and locally.

**Trade-off:** docs are served only when the CoolSolve binary is running; they are not independently accessible on the internet without the dedicated server also running CoolSolve.

---

## Option E — Self-hosted dedicated server: Docusaurus / VitePress

A full-featured documentation framework (e.g. [VitePress](https://vitepress.dev/) or [Docusaurus](https://docusaurus.io/)) hosted on the dedicated server, built from the `.md` files by a CI pipeline (GitHub Actions pushes the built site to the server via SSH or `rsync`).

**Dependencies introduced:**
- Node.js + npm (only in CI)
- VitePress or Docusaurus (npm packages, locked to a version in `package.json`)
- A GitHub Actions deploy key to the dedicated server
- Node must be kept up to date as the framework evolves

**Advantages over A:** richer UI, versioning, algolia search, custom domain, full control.

**Trade-off:** heaviest maintenance burden; npm ecosystem churn; requires keeping the dedicated server running a Node.js-built static site or a reverse proxy.

---

## Option F — Hybrid: GitHub Pages (Option A or C) + in-GUI deep link

Any of the above options for GitHub Pages is combined with a single "Open documentation" link in the CoolSolve GUI toolbar that points to the GitHub Pages URL. The dedicated server can optionally mirror or redirect there. No duplication, no extra code in the binary.

This is not a standalone option but a composition layer applicable on top of A, B, or C.

---

## Summary Table

| Option | Source of truth | Build step required | External service | Offline docs | Search | Custom nav | Maintenance burden | Dependencies |
|--------|----------------|--------------------|-----------------|-----------|---------|-----------|--------------------|--------------|
| **A — MkDocs + GH Pages** | `.md` files | Yes (Python, CI) | GitHub Pages | No | Yes (lunr.js built-in) | Yes (mkdocs.yml) | Low–Medium | `mkdocs`, GitHub Actions |
| **B — Pandoc + GH Pages** | `.md` files | Yes (shell, CI) | GitHub Pages | No | No | Minimal (hand-written) | Low | `pandoc`, GitHub Actions |
| **C — Jekyll GH Pages** | `.md` files | None (GitHub-side) | GitHub Pages | No | No | Minimal | Very low | None |
| **D — Embedded in binary** | `.md` files | Yes (C++ embed at build) | None | Yes | No | Minimal (JS) | Low | `marked.js` (one file) |
| **E — Dedicated server (VitePress)** | `.md` files | Yes (Node, CI+deploy) | Dedicated server | No | Yes (Algolia/built-in) | Yes | High | Node, npm, SSH deploy |
| **F — GH Pages + GUI link** | `.md` files | Depends on A/B/C | GitHub Pages | Partial (D hybrid) | Depends | Depends | Low | Depends on base option |

---

## Recommended Solution: **Option A (MkDocs-Material) + Option F (GUI deep link)**

### Rationale

**Single source of truth is already satisfied** — the `.md` files are never duplicated; MkDocs reads them directly.

**MkDocs is the lightest framework that still provides:**
- A navigation sidebar defined once in `mkdocs.yml`
- Full-text search (built-in, client-side, no Algolia account)
- Syntax highlighting for code blocks
- A professional appearance (MkDocs-Material theme)
- `pip install mkdocs-material` — one command, no npm, no Node.js
- The Python dependency only runs in CI; it is **never part of the C++ build**, so it does not affect the binary or its dependencies

**The GitHub Actions workflow** is ~15 lines and runs only on `main` pushes. It requires no secrets and no external accounts beyond a standard GitHub repo.

**The GUI link** (Option F) is a single `<a href="...">` in the toolbar — three lines of TypeScript, no new C++ code, no new dependencies. The URL is the permanent GitHub Pages address, which works even when the binary is not running.

**The dedicated server**, if desired later, can simply serve the same pre-built static site (rsync from the `gh-pages` branch) with zero additional tooling.

### What to add to the repository

```
mkdocs.yml                    ← navigation and theme config (new, ~20 lines)
.github/
  workflows/
    deploy-docs.yml           ← CI workflow to build and deploy (new, ~15 lines)
docs/
  syntax_reference.md         ← new page: CoolSolve/EES language syntax
  index.md                    ← new: landing page / table of contents
  (all existing .md files unchanged)
gui/src/components/Toolbar.tsx  ← add "Documentation" link (3 lines)
```

### What is explicitly avoided

- No HTML files to maintain alongside the `.md` files
- No Node.js / npm in the documentation pipeline
- No Docusaurus, no VitePress, no Jekyll configuration beyond a `_config.yml` fallback
- `marked.js` (Option D) is not needed because the GUI will simply link out to GitHub Pages — the binary stays lean

### Upgrading later

If the documentation outgrows GitHub Pages (e.g. needs versioning or a custom domain on the dedicated server), the switch from GitHub Pages + MkDocs to a self-hosted MkDocs deployment requires changing exactly one line in the CI workflow (`mkdocs gh-deploy` → `rsync dist/ user@server:/var/www/docs`). No content changes needed.
