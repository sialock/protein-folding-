# Contributing

We welcome reproducibility, correctness, test coverage and documentation
improvements. Before submitting a pull request:

```bash
python -m pip install -e ".[dev]"
python -m pytest -q
python -m compileall -q src
```

- Add tests for numerical changes and validation failures.
- Specify physical units for all future energy terms and simulation parameters.
- Do not state that a reference-biased toy calculation predicts protein folding.
- Keep simulation inputs local and tests independent of network access.
- Prefer narrowly scoped pull requests documenting the expected scientific effect.