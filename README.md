# colocalized-genes-chart

## requirements

The Data-Driven Ontology Toolkit (DDOT) - https://github.com/idekerlab/ddot
NetColoc - https://github.com/ucsd-ccbb/NetColoc

## setup

DDOT depends on Python libs that need Python <=3.9

```
git clone https://github.com/Fuzzwah/colocalized-genes-chart.git
cd colocalized-genes-chart
python3.9 -m venv .env
. .env/bin/activate
pip install -U pip
pip install -r requirements.txt
git clone --branch python3 https://github.com/idekerlab/ddot.git
cd ddot
python setup.py build
python setup.py install
cd ..
```

## running

Drop your list of genes into the `input` directory, then run `python main.py`. You'll be prompted with a list of the files in the `input` directory, select the one you'd like to process