# Design plans (not served)

This folder holds design plans and their prototypes, for contributors and coding agents. It is a
subfolder of `docs/`, and `cmake/embed_assets.cmake` embeds only the files directly in `docs/`
(subfolders are skipped), so nothing here is served by a CoolSolve server under `/docs/`. Do not
link these files from the served documentation (`README.md`, `docs/*.md`): the links would be
broken in the in-app viewer.

| Plan | Status |
|---|---|
| [Dashboards](dashboards.md) — HTML dashboards of a model: values on schematics, inputs, GUI buttons, saved diagrams as tabs, animation of `INTEGRAL` models, LLM authoring; prototype and examples in [`dashboards/`](dashboards/) | proposal, 2026-11-20 |
