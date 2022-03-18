# How to build the documentation with sphinx

1. From this directory `docs-src` on your local machine, do 

```make html ```

This will build the html documentation into the directory `docs-src/build/html`. Note that the `docs-src/build` directory is **not** sync to the github repository (i.e. it is in the .gitignore). This is so that you can build the documentation locally and check it before you commit it and update the web-hosted documentation. 

2. To check your changes, you can locally open `docs-src/build/html/index.html` to check the html version of the documentation.

2. If you are satisfied with your changes, then you do 

```make github ```

This does a `make html` and then copy the content of the `docs-src/build/html` into `docs`, which is the folder in which the github pages are hosted. The online-version of the documentation should then be udated after you push and commit your modifications. 

## Local pre-requisites to build the documentation locally.
