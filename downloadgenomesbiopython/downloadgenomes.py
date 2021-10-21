import urllib

import tensorflow_decision_forests as tfdf

import os
import numpy as np
import pandas as pd

try:
    from wurlitzer import sys_pipes
except Exception as exp:
    import colabtools


def split_dataset(dataset, test_ratio=0.30):
    """Splits a panda dataframe in two."""
    test_indices = np.random.rand(len(dataset)) < test_ratio
    return dataset[~test_indices], dataset[test_indices]


def findgenesofsignificance():
    # !wget - q https: // storage.googleapis.com / download.tensorflow.org / data / palmer_penguins / penguins.csv - O / tmp / penguins.csv

    # Load a dataset into a Pandas Dataframe.
    dataset_df = pd.read_csv("/data/MicroarrayExpression.csv")

    # Display the first 3 examples.
    dataset_df.head(3)

    train_ds_pd, test_ds_pd = split_dataset(dataset_df)
    print("{} examples in training, {} examples for testing.".format(
        len(train_ds_pd), len(test_ds_pd)))

    train_ds = tfdf.keras.pd_dataframe_to_tf_dataset(train_ds_pd)
    test_ds = tfdf.keras.pd_dataframe_to_tf_dataset(test_ds_pd)

    # Specify the model.
    model_1 = tfdf.keras.RandomForestModel()

    # Optionally, add evaluation metrics.
    model_1.compile(
        metrics=["accuracy"])

    # Train the model.
    # "sys_pipes" is optional. It enables the display of the training logs.
    with sys_pipes():
        model_1.fit(x=train_ds)
    evaluation = model_1.evaluate(test_ds, return_dict=True)
    print()

    for name, value in evaluation.items():
        print(f"{name}: {value:.4f}")

    tfdf.model_plotter.plot_model_in_colab(model_1, tree_idx=0, max_depth=3)
    model_1.summary()
    model_1.make_inspector().features()
    model_1.make_inspector().variable_importances()
    model_1.make_inspector().evaluation()
    model_1.make_inspector().training_logs()
    model_1.make_inspector().export_to_tensorboard("/tmp/tensorboard_logs")


def get_assembly_summary(identifier: str):
    """Get esummary for an entrez id"""
    from Bio import Entrez
    esummary_handle = Entrez.esummary(db="assembly", id=identifier, report="full")
    esummary_record = Entrez.read(esummary_handle)
    return esummary_record


def get_assemblies(term, download=True, path='assemblies'):
    """Download genbank assemblies for a given search term.
    Args:
        term: search term, usually organism name
        download: whether to download the results
        path: folder to save to
    """

    from Bio import Entrez
    # provide your own mail here
    Entrez.email = "lalitha.viswanathan79@gmail.com"
    handle = Entrez.esearch(db="assembly", term=term, retmax='200')
    record = Entrez.read(handle)
    ids = record['IdList']
    print(f'found {len(ids)} ids')
    links = []
    for ident in ids:
        # get summary
        summary = get_assembly_summary(ident)
        # get ftp link
        url = summary['DocumentSummarySet']['DocumentSummary'][0]['FtpPath_RefSeq']
        if url == '':
            continue
        label = os.path.basename(url)
        link = os.path.join(url, label + '_genomic.fna.gz')
        print(link)
        links.append(link)
        if download:
            urllib.request.urlretrieve(link, f'{label}.fna.gz')

    return links
