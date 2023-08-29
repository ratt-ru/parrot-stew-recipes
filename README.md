# parrot-stew-recipes

Recipes etc. for the RATT PARROT data reduction 

## 1GC 

You will need to install CARACal 1.0.6+ (https://github.com/caracal-pipeline/caracal) for the 1GC processing -- see ``caracal-configs`` for config files.

Once 1GC has been done and the target MS has been extracted, proceed to the selfcal/lightcurve recipes

## Selfcal, dynamic imaging and/or lightcurve extraction

The entire post-1GC workflow is implemented as a Stimela recipe (https://stimela.readthedocs.io), and should be fully reproducible given a big enough compute node with a Singularity installation.

Prerequisites:

*   stimela 2.0 (https://github.com/caracal-pipeline/stimela), tag ``parrot1``, or version 2.0rc8

*   cult-cargo (https://github.com/caracal-pipeline/cult-cargo), tag ``parrot1``, or version 0.2.0

    Since stimela 2 and the associated cult-cargo package were in pre-release at time of publication, we specify a git tag here rather than a version. Subsequent versions may also work, but are not guaranteed to, given the code churn inherent to pre-releases.

    To install the tagged versions, create a virtual environment, then do e.g.:

    ```
    pip install git+https://github.com/caracal-pipeline/stimela.git@parrot1
    pip install git+https://github.com/caracal-pipeline/cult-cargo.git@parrot1
    ```

*   The ``omstimelation`` (https://github.com/o-smirnov/omstimelation) recipe collection, tag ``parrot1``

*   This recipe collection (https://github.com/ratt-ru/parrot-stew-recipes), tag ``parrot1``.

    The collections aren't formal packages (at time of writing) and don't need to be installed, just checked out from github:

    ```
    git clone https://github.com/ratt-ru/parrot-stew-recipes -b parrot1
    cd parrot-stew-recipes
    git clone https://github.com/o-smirnov/omstimelation -b parrot1
    ```






