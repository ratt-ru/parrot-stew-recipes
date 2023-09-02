# parrot-stew-recipes

Recipes etc. for the RATT PARROT data reduction 

## 1GC 

You will need to install CARACal 1.0.6+ (https://github.com/caracal-pipeline/caracal) for the 1GC processing -- see ``caracal-configs`` for config files.

Once 1GC has been done and the target MS has been extracted, proceed to the selfcal/lightcurve recipes

## Recipe installation

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

    The collections are not formal packages (at time of writing) and don't need to be (and can't be) pip installed. Rather, just check out appropriate tags/branches from github as follows:

    ```
    git clone https://github.com/ratt-ru/parrot-stew-recipes -b parrot1
    cd parrot-stew-recipes
    git clone https://github.com/o-smirnov/omstimelation -b parrot1
    ```

You're now all set to go.

## Selfcal & dynamic imaging 

The conjunction observation needs to be processed scan by scan. There is a "preparation" recipe, which splits out one scan into a per-scan MS, and updates field/UVW coordinates:

```
$ stimela run jove-prepare.yml scan=4
```

Then there is the a selfcal recipe, which does imaging and selfcal on one (previously prepared) per-scan MS:

```
$ stimela run jove-pol.yml scan=4
```

Finally, there is a wrapper recipe, which can be used to execute the above two recipes in a loop over multiple scans:

```
$ stimela run jove-pol-loop.yml                                       # prepares and selfcals all scans 
$ stimela run jove-pol-loop.yml scan-list=[4, 6, 8]                   # prepares and selfcals specific scans
$ stimela run jove-pol-loop.yml scan-list=[4, 6, 8] -s jove-prepare   # prepares speciic scans
$ stimela run jove-pol-loop.yml scan-list=[4, 6, 8] -s jove-pol       # runs selfcal on speciic scans
```

## Selfcal and lightcurve extraction for the follow-up observations

The recipe uses static deconvolution masks. Extract them first using ``tar zxvf masks.tgz``.

Run the imaging recipe, specifying a particular observation via the ``obs`` input:

```
$ stimela run image-parrot.yml obs=U2
```

See ``parrot-observation-sets.yml`` for a list of observations.

A copy of ``mastercat.ecsv`` is provided in the repository. This catalog can also be created via the ``make-master-catalog`` step, based of source finder outputs from multiple observations.


