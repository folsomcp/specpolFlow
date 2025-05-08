# Contributing

SpecpolFlow is open source software, developed by a community of volunteers.  And we could use your help!  

There are several ways you can contribute to the project:
* Report bugs, or potential errors -- We try to test the things carefully, but in any sufficiently complicated code mistakes can slip through, so let us know what problems you run into.
* Suggestions -- If there is a feature that you would like to have, or a way to make SpecpolFlow easier to use, let us know.  We can't promise to implement every suggestion (we have limited time!) but we'd at least like to know what would be useful.
* Contribute code -- If there is a feature you have, this is a good way to make sure it happens!  
* Add documentation -- A software package is only as good as its documentation.  We try to go beyond the bare minimum, and show different ways to use SpecpolFlow, for people with different levels of experience.  So there is always room to improve or expand the documentation and tutorials.

## Reporting issues and making suggestions

If you report a problem, or if you just want to make a suggestion, an easy way to do that is through the GitHub [issue tracker](https://github.com/folsomcp/specpolFlow/issues).  You can create a new thread (with the 'New issue' button), and also see what issues people have been talking about.

You can also contact us directly at specpolflow@gmail.com

## Contributing code

If you want to get serious about contributing code, then the a good place to start is cloning the GitHub repository.  However, if you're not familiar with Git, and just want to contribute a small bit of code or a stand alone function, you can post it in the GitHub [Issues page](https://github.com/folsomcp/specpolFlow/issues), or email us at specpolflow@gmail.com

Inside the Git repository, the documentation is all in the `docs-jb` folder, with folders for the `About`, `GetStarted`, and `Tutorials` sections of the website.  The documentation is all in the form of Markdown files, and Jupyter notebooks.  The Python code for the project is all in the `specpolFlow` folder.

### Setting up with Git

To make a copy of SpecpolFlow that you can modify, the typical way to do that is to start by [forking](https://docs.github.com/en/pull-requests/collaborating-with-pull-requests/working-with-forks/about-forks) the repository.  This basically creates your own copy of the code on GitHub.  You can do that on the [GitHub page](https://github.com/folsomcp/specpolFlow) with the `fork` button at the top right.  

You can then make a local copy of your version of SpecpolFlow using the git clone command, for example:
```
git clone https://github.com/YOUR-GITHUB-USERNAME/specpolFlow.git
```
That will download the code and documentation files, which you can look at and modify.

If you want to install your local copy of SpecpolFlow, it is a good idea to first make a [virtual environment](https://docs.python.org/3/tutorial/venv.html) to install the code in, using either Python's venv or conda.  For example create a virtual environment called `spf-dev` with venv:
```
python -m venv spf-dev-env
```
That will create a new folder `spf-dev-env`.  Then activate the virtual environment by running the activate script:
```
source spf-dev-env/bin/activate
```
This keeps your development version separate from the regular version of SpecpolFlow.  You can switch versions just by activating or deactivating (just type `deactivate`) this virtual environment.  If you are using an IDE like VS Code, it may be easier to set up this virtual environment from within the that program.  

You can make an ['editable' install](https://setuptools.pypa.io/en/latest/userguide/development_mode.html) of your copy of SpecpolFlow with pip.  First go to the folder that has your copy of SpecpolFlow (it should have pyproject.toml and README.md), then run pip with the `-e` flag:
```
cd specpolFlow
pip install -e .
```
(The . at the end tells pip to install the version 'here' in that specpolFlow folder, not from their servers.)

Now you can edit the code and test the modified code without having to uninstall and reinstall every time you change something.  Once you have a set of changes you are happy with, you can update Git with them.  First you need to `add` the file that you modified or created.  Then you need to `commit` the modified files, usually including a note about what you have changed.  
```
git add specpolFlow/file_that_has_been_modified.py
git commit -m "a helpful message describing the changes you made"
```
That saves the changes in your local copy of Git.  Then you should update your version on GitHub to match the local version using `push`.
```
git push 
```
For more information about Git, GitHub has some [beginner](https://docs.github.com/en/get-started/git-basics) and [more extensive](https://docs.github.com/en/get-started/using-git), and the Git project has extensive [documentation](https://git-scm.com/docs).

You can make multiple commits, and push changes to your fork on GitHub multiple times, as you develop and test.  If you change more than one thing, it is often better to make those changes as separate commits.  If you change how something functions, or add new functions, please also update the documentation!

Once you are done testing, and have made all the changes you want to make, you can [pull request](https://docs.github.com/en/pull-requests/collaborating-with-pull-requests/proposing-changes-to-your-work-with-pull-requests/creating-a-pull-request-from-a-fork) using the GitHub website.  By making a pull request, you are asking to take the changes from your fork (your copy) of the code and add them to the main version of SpecpolFlow.  

### Building a local copy of the website

If you want to create a local copy of the website, you first need a local copy of all the code, so follow the above Setting up with Git section.  

The website is built with the [Jupyter Book](https://jupyterbook.org/) package.  If you don't yet have it, Jupyter Book can be installed with pip or conda, e.g.:
```
pip install jupyter-book
```

You can then run Jupyter Book from the main SpecpolFlow folder to [build the book](https://jupyterbook.org/en/stable/start/build.html).  You should see the `docs-jb` folder there, and you can point the `jupyter-book build` command at that folder:
```
jupyter-book build docs-jb/
```

This will generate a folder inside `docs-jb` called `_build`, with the website in `docs-jb/_build/html/`, and the front page of the website should be `docs-jb/_build/html/index.html`.  
