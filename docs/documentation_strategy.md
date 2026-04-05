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
- **Syntax reference** — a new page documenting the CoolSolve language must be added

---

Different options were considered:
- Option A — GitHub Pages with MkDocs (static site, CI-generated)
- 

- Option B — GitHub Pages with Pandoc (no Python framework)
-Option C — GitHub Pages with no build step (plain Markdown rendering)
- Option D — Self-hosted: serve the `.md` files via the CoolSolve binary itself
- Option E — Self-hosted dedicated server: Docusaurus / VitePress
- Option F — Hybrid: GitHub Pages (Option A or C) + in-GUI deep link
## Selected Option: Self-hosted: serve the `.md` files via the CoolSolve binary itself

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


### What is explicitly avoided

- No HTML files to maintain alongside the `.md` files
- No Node.js / npm in the documentation pipeline
- No Docusaurus, no VitePress, no Jekyll configuration beyond a `_config.yml` fallback
- `marked.js` (Option D) is not needed because the GUI will simply link out to GitHub Pages — the binary stays lean

### Upgrading later

If the documentation outgrows GitHub Pages (e.g. needs versioning or a custom domain on the dedicated server), the switch from GitHub Pages + MkDocs to a self-hosted MkDocs deployment requires changing exactly one line in the CI workflow (`mkdocs gh-deploy` → `rsync dist/ user@server:/var/www/docs`). No content changes needed.
