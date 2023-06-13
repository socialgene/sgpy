# Contributing

Contributions are welcome, and they are greatly appreciated! Every
little bit helps, and credit will always be given.

## Bug reports

When [reporting a bug](https://chasemc/socialgene/bin/socialgene/issues) please include:

    * Your operating system name and version.
    * Any details about your local setup that might be helpful in troubleshooting.
    * Detailed steps to reproduce the bug.

## Documentation improvements

SocialGene could always use more documentation, whether as part of the [official docs](https://github.com/socialgene/socialgene.github.io), or  even on the web in blog posts, articles, etc.

## Feature requests and feedback

The best way to send feedback is to file an issue at <https://github.com/socialgene/sgpy/issues/new/choose>

If you are proposing a feature:

* Explain in detail how it would work.
* Keep the scope as narrow as possible, to make it easier to implement.
* Remember that this is a volunteer-driven project, and that code contributions are welcome :)

## Development

To set up this Python package for local development (requires git, Make, and conda (optional)):

1. Navigate to [socialgene/sgpy](https://github.com/socialgene/sgpy) and "Fork" the repo (look for the "Fork" button).
2. Clone your fork locally:

```bash
      git clone git@github.com:socialgene/sgpy.git
```

3. Create a branch for local development::

```bash
    git checkout -b name-of-your-bugfix-or-feature
```

   Now you can make your changes locally.

4. When you're done making changes run all the checks using make

```bash
    make run_ci
```

5. Commit your changes and push your branch to GitHub::
  <https://www.conventionalcommits.org/en/v1.0.0-beta.2/>

    git add .
    git commit -m "Your detailed description of your changes."
    git push origin name-of-your-bugfix-or-feature

6. Submit a pull request through the GitHub website.

## Pull Request Guidelines

If you need some code review or feedback while you're developing the code just make the pull request.

For merging, you should:

1. Include passing tests (run ``tox``) [1]_.
2. Update documentation when there's new API, functionality etc.
3. Add a note to ``CHANGELOG.rst`` about the changes.
4. Add yourself to ``AUTHORS.rst``.

## Tips

To run a subset of tests::

    tox -e envname -- pytest -k test_myfeature

To run all the test environments in *parallel*::

    tox -p auto
