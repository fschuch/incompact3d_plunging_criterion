# Building these Docs

This website is built with [Jupyter Book](https://jupyterbook.org)
from the [`docs`](https://github.com/fschuch/incompact3d_plunging_criterion/tree/main/docs)
folder of the repository. The dependencies are declared in
[`docs/pyproject.toml`](https://github.com/fschuch/incompact3d_plunging_criterion/blob/main/docs/pyproject.toml)
and pinned by
[`docs/uv.lock`](https://github.com/fschuch/incompact3d_plunging_criterion/blob/main/docs/uv.lock)
for reproducible builds; [uv](https://docs.astral.sh/uv/) manages the
environment automatically in `docs/.venv`.

With uv installed, build from the repository root:

```bash
make docs        # one-off build into docs/_build/html
```

```bash
make docs-serve  # live-reloading preview in the browser
```

```{important}
The example notebooks are rendered with their stored outputs, since
the dataset they depend on is not included in the repository (see the
[dataset description](../dataset/description.md)). They are never
executed at build time, so no scientific software stack is required
to build the documentation.
```

## Continuous deployment

The GitHub Actions workflow
[docs.yaml](https://github.com/fschuch/incompact3d_plunging_criterion/blob/main/.github/workflows/docs.yaml)
builds the book with warnings treated as errors on every pull request
touching the documentation, and deploys it to GitHub Pages on every
push to the `main` branch.
