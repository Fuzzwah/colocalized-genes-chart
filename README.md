# colocalized-genes-chart

## requirements

- The Data-Driven Ontology Toolkit (DDOT) - https://github.com/Fuzzwah/ddot.git
- NetColoc - https://github.com/ucsd-ccbb/NetColoc

## setup

Note: DDOT depends on Python libs that need Python <=3.8

```
git clone https://github.com/Fuzzwah/colocalized-genes-chart.git
cd colocalized-genes-chart
python3.8 -m venv .env
. .env/bin/activate
pip install -U pip
pip install -r requirements.txt
git clone https://github.com/Fuzzwah/ddot.git
cd ddot
python setup.py build
python setup.py install
cd ..
```

## running

Drop your list of genes into the `input` directory, then run `python main.py`. You'll be prompted with a list of the files in the `input` directory, select the one you'd like to process.

## example output

![image](https://github.com/Fuzzwah/colocalized-genes-chart/assets/658935/36f9947a-b493-47f9-a6d2-b82ee6679b6a)
