# C2VSIM-FG Stream Transects

Stream transects created using C2VSIM-FG, CASGEM, and AEM lithology logs. This repository uses Poetry for version control

### Building the book

If you'd like to develop and/or build the C2VSIM-FG Stream Transects book, you should:

1. Clone this repository
2. Run `poetry install` within the directory containing the repository
3. (Optional) Edit the books source files located in the `c2vsim_fg_stream_transects/` directory
4. Since we only track the .py files, you need to build the .ipynb twins using jupytext. Run `poetry run jupytext --set-formats ipynb,py:percent c2vsim_fg_stream_transects/*.py` when you clone the repository, and
   `poetry run jupytext --update c2vsim_fg_stream_transects/*.py` when you make edits.
6. Run `poetry run jupyter-book clean c2vsim_fg_stream_transects/` to remove any existing builds
7. Run `poetry run jupyter-book build c2vsim_fg_stream_transects/`

A fully-rendered HTML version of the book will be built in `c2vsim_fg_stream_transects/_build/html/`.

### Hosting the book

Please see the [Jupyter Book documentation](https://jupyterbook.org/publish/web.html) to discover options for deploying a book online using services such as GitHub, GitLab, or Netlify.

For GitHub and GitLab deployment specifically, the [cookiecutter-jupyter-book](https://github.com/executablebooks/cookiecutter-jupyter-book) includes templates for, and information about, optional continuous integration (CI) workflow files to help easily and automatically deploy books online with GitHub or GitLab. For example, if you chose `github` for the `include_ci` cookiecutter option, your book template was created with a GitHub actions workflow file that, once pushed to GitHub, automatically renders and pushes your book to the `gh-pages` branch of your repo and hosts it on GitHub Pages when a push or pull request is made to the main branch.

## Contributors

We welcome and recognize all contributions. You can see a list of current contributors in the [contributors tab](https://github.com/jtdiazcLWA/c2vsim_fg_stream_transects/graphs/contributors).

## Credits

This project is created using the excellent open source [Jupyter Book project](https://jupyterbook.org/) and the [executablebooks/cookiecutter-jupyter-book template](https://github.com/executablebooks/cookiecutter-jupyter-book).
