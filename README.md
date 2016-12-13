# nextHIV

## Description

Processing code for nextHIV. The package skeleton was created using [makenew/python-package](https://github.com/makenew/python-package).

## Docker installation

## Local installation

### Source Code

The [nexthiv source](https://github.com/sdwfrost/nexthiv) is hosted on
GitHub. Clone the project with

```bash
git clone https://github.com/sdwfrost/nexthiv.git
```

### Requirements

You will need [Python 3](https://www.python.org/) with
[pip](https://pip.pypa.io/).

Install the dependencies with

```bash
pip install -r requirements.txt
```

## Usage

### Configuration

Configuration for the database is stored as a YAML document, An example is found below.

```yaml
rethinkdb:
  host: 'localhost'
  port: 28015
  db: 'nexthiv'
tables:
  - sequences
  - processed_sequences
  - demographics
  - regimens
  - labs
  - references
  - resistance
  - subtypes
  - clustering
  - phylogenies
sequence:
  table: 'sequences'
  name: 'SEQUENCE'
  startpos: 0
  endpos: 1497
  processed_table: 'processed_sequences'
  processed_name: 'SEQUENCE_ALIGNED'
directories:
  datadir: './data'
  tmpdir: './tmp'
programs:
  mafft:
    cmd: 'mafft'
  tn93:
    cmd: 'tn93'
    threshold: 1.0
    ref_threshold: 0.05
    ambiguities: 'resolve'
    min_overlap: 100
    fraction: 1.0
  iqtree:
    cmd: 'iqtree-omp'
    threads: 4
    model: 'GTR+G4'
    bootstrap_samples: 1000
```

### Initialising database

```bash
nexthiv init -c nexthiv.yaml
```

### Importing data

```bash
nexthiv import -c nexthiv.yaml -i sequences.txt -d sequences
```

### Aligning sequences

The following command will align the sequences against a reference, trim them, and insert them into the `processed_sequences` table.

```bash
nexthiv align -c nexthiv.yaml
```

### Clustering sequences

```bash
nexthiv cluster -c nexthiv.yaml
```

## Contributing

Please submit and comment on bug reports and feature requests.

To submit a patch:

1.  Fork it (<https://github.com/sdwfrost/nexthiv/fork>).
2.  Create your feature branch (`git checkout -b my-new-feature`).
3.  Make changes. Write and run tests.
4.  Commit your changes (`git commit -am 'Add some feature'`).
5.  Push to the branch (`git push origin my-new-feature`).
6.  Create a new Pull Request.

## License

This Python package and associated files are licensed under the MIT license.

## Warranty

This software is provided "as is" and without any express or implied
warranties, including, without limitation, the implied warranties of
merchantibility and fitness for a particular purpose.
