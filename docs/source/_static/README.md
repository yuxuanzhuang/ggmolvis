# Static Doc Directory

Add any paths that contain custom static files (such as style sheets) here,
relative to the `conf.py` file's directory. 
They are copied after the builtin static files,
so a file named "default.css" will overwrite the builtin "default.css".

The path to this folder is set in the Sphinx `conf.py` file in the line: 
```python
html_static_path = ['_static']
```

## Examples of file to add to this directory
* Custom Cascading Style Sheets
* Custom JavaScript code
* Static logo images


### Logos

Several template and placeholder logo documents are
already provided in this repo.
We encourage you to replace them with your own custom logos!
You can either use an entirely new image, or edit our
template documents (for example, by placing your own logo
within the gears).

All "Empty gear" and placeholder documents in the ``logo/`` directory are provided
under a [CC-BY-4.0](https://creativecommons.org/licenses/by/4.0/) license under the [terms of our license document](https://github.com/MDAnalysis/branding/blob/main/logos/LICENSE).
