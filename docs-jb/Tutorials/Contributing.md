# Contributing and reporting issues

SpecpolFlow is open source software developed by a community of users. We could use your help!  

There are several ways you can contribute to the project:
* Report bugs, or potential errors -- We try to test things carefully, but mistakes can slip through, so let us know if you run into problems.
* Suggestions -- If there is a feature that you would like to have, or a way to make SpecpolFlow easier to use, let us know.  We can't promise to implement every suggestion (we have limited time!), but we'd like to know what would be useful.
* Contribute code -- If there is a feature you'd like to have, this is a good way to make sure it happens!  If you find a solution to a bug, please share it with the community.
* Add documentation -- A software package is only as good as its documentation.  We try to go beyond the bare minimum, with tutorials on using SpecpolFlow for people with different levels of experience.  There is always room to improve or expand the documentation and tutorials.

If you would like to contribute to SpecpolFlow on a regular basis, contact us about meeting with the core development team.

## Reporting issues and making suggestions

If you want to report a problem, or just make a suggestion, an easy way to do that is through the GitHub [issue tracker](https://github.com/folsomcp/specpolFlow/issues).  You can create a new thread (with the 'New issue' button), and also see what issues people have been talking about.

You can also contact us directly at specpolflow@gmail.com.

## Contributing documentation

The documentation for SpecpolFlow is all in the form of Jupyter notebooks and Markdown files.  These are all in the [`docs-jb`](https://github.com/folsomcp/specpolFlow/tree/main/docs-jb) folder on GitHub.  Inside that are folders for the `About`, `GetStarted`, and `Tutorials` sections of the website.

If you want to make smaller changes, you can just download a copy of the relevant file, edit it, and contact us through email or on the [Issues page](https://github.com/folsomcp/specpolFlow/issues).  If you want to make bigger changes across multiple files, or are just comfortable using Git, you can follow the directions in the 'Contributing code' section below to download all of the code and documentation files.

:::{note}
The .md and Jupyter files contain some non-standard Markdown notation, using extensions from [MyST Markdown](https://myst-parser.readthedocs.io/).  This is used for rendering some features on the website with [Jupyter Book](https://jupyterbook.org/) (like this [note](https://jupyterbook.org/en/stable/content/content-blocks.html#notes-warnings-and-other-admonitions)!).
:::

## Contributing code

If you want to get serious about contributing code, a good place to start is forking the GitHub repository.  However, if you're not familiar with Git and just want to contribute a small bit of code or a stand alone function, you can post it in the GitHub [Issues page](https://github.com/folsomcp/specpolFlow/issues), or email us at specpolflow@gmail.com.

The Python code for the project is all in the `specpolFlow` folder, while documentation and tutorials are in the `docs-jb` folder.

## Setting up with Git

### Fork the GitHub repository

To make a copy of SpecpolFlow that you can modify, the typical way is to start with [forking](https://docs.github.com/en/pull-requests/collaborating-with-pull-requests/working-with-forks/about-forks) the repository.  This basically creates your own copy of the code on GitHub.  You can do that on the SpecpolFlow [GitHub page](https://github.com/folsomcp/specpolFlow) with the `fork` button at the top right.  

You can then download a local copy of your version of SpecpolFlow using the git clone command, for example:
```
git clone https://github.com/YOUR-GITHUB-USERNAME/specpolFlow.git
```
That will download the code and documentation files, which you can look at and modify.

Optionally, you may wish to set up your fork so that it can be updated if there are any changes to the upstream repository (the main version of SpecpolFlow).  First [define the upstream repository](https://docs.github.com/en/pull-requests/collaborating-with-pull-requests/working-with-forks/configuring-a-remote-repository-for-a-fork), then download any changes [synching your fork](https://docs.github.com/en/pull-requests/collaborating-with-pull-requests/working-with-forks/syncing-a-fork#syncing-a-fork-branch-from-the-command-line).

### Make a virtual environment

If you want to install your local copy of SpecpolFlow, it is a good idea to first make a [virtual environment](https://docs.python.org/3/tutorial/venv.html) to install the code in, using either venv or conda.  The idea is that packages installed in one virtual environment are only used inside that virtual environment.  This lets you keep your development version of SpecpolFlow separate from the regular version, and keeps the development version out of your regular Python environment.  For example, to create a virtual environment called `spf-dev-env` with venv:
```
python -m venv spf-dev-env
```
That will create a new folder `spf-dev-env`.  Then activate the virtual environment by running a script, on Mac and Linux run:
```
source spf-dev-env/bin/activate
```
or on Windows run:
```
spf-dev-env\Scripts\activate
```
With the virtual environment active, you can now install the development version of SpecpolFlow into this environment, instead of your regular Python environment.  If you have a different version installed in different environments, you can switch versions just by activating or deactivating (just type `deactivate`) the virtual environment.  Note, if you are using an integrated development environment (IDE) like VS Code, it may be easier to set up a virtual environment from within that program.  

### Install your local copy

You can make an [editable install](https://setuptools.pypa.io/en/latest/userguide/development_mode.html) of your copy of SpecpolFlow with pip.  This installs the package while still letting you edit the source files in their current location.  First go into the folder that has your copy of SpecpolFlow (it should have pyproject.toml and README.md), then run pip with the `-e` flag:
```
cd specpolFlow
pip install -e .
```
(The `.` at the end tells pip to install the version 'here' in that specpolFlow folder.)

### Save your changes with Git

Now you can edit and test the code without having to uninstall and reinstall every time you change something.  Once you have a set of changes you are happy with, you can update Git with them.  First you need to `add` the file that you modified or created.  Then you need to `commit` the modified files, usually including a note about what you have changed.  
```
git add specpolFlow/file_that_has_been_modified.py
git commit -m "a helpful message describing the changes you made"
```
That saves the changes to your local Git archive.  Then you should update your version on GitHub to match the local version using `push`.
```
git push 
```
For more information about Git, GitHub has some [beginner](https://docs.github.com/en/get-started/git-basics) and [more extensive](https://docs.github.com/en/get-started/using-git) documentation, and the Git project has [extensive](https://git-scm.com/docs) documentation.

You can make multiple commits, and push the changes to your fork on GitHub, as you develop and test your work.  Often, if you change more than one thing, it is better to make those changes as separate commits.  Also, when you change how something functions (or add new functions) please update the documentation!

### Submit a pull request

Once you are done testing, and have made all the changes you want to make, you can submit a [pull request](https://docs.github.com/en/pull-requests/collaborating-with-pull-requests/proposing-changes-to-your-work-with-pull-requests/creating-a-pull-request-from-a-fork) using the GitHub website.  By making a pull request, you are asking to take the changes from your fork (your copy) of the code and add them to the main version of SpecpolFlow.  


## Building a local copy of the website

If you want to create a local copy of the website, you first need a local copy of all the code, so see the above 'Setting up with Git' section.  

The website is built with the [Jupyter Book](https://jupyterbook.org/) package.  If you don't yet have it, Jupyter Book can be installed with pip or conda, e.g.:
```
pip install jupyter-book
```

You can then run Jupyter Book from the main SpecpolFlow folder to [build the book](https://jupyterbook.org/en/stable/start/build.html).  You should see the `docs-jb` folder there, and you can point the `jupyter-book build` command at that folder:
```
jupyter-book build docs-jb/
```

This will generate a folder inside `docs-jb` called `_build`, with the website in `docs-jb/_build/html/`, and the front page of the website should be `docs-jb/_build/html/index.html`.  
