from setuptools import setup, find_packages

setup(
    name='cortado_core',
    version='1.3.6',
    author="Daniel Schuster",
    author_email="daniel.schuster@fit.fraunhofer.de",
    # package_dir={"": "cortado-core"},
    # packages=find_packages(where="cortado_core"),
    packages=['cortado_core', 'cortado_core.utils','cortado_core.models', 'cortado_core.freezing', 'cortado_core.process_tree_utils',
              'cortado_core.performance', 'cortado_core.variant_query_language',
              'cortado_core.variant_query_language.grammars', 'cortado_core.alignments',
              'cortado_core.alignments.infix_alignments', 'cortado_core.alignments.infix_alignments.variants',
              'cortado_core.alignments.prefix_alignments', 'cortado_core.alignments.prefix_alignments.variants',
              'cortado_core.alignments.suffix_alignments','cortado_core.subprocess_discovery', 'cortado_core.subprocess_discovery.subtree_mining',
              'cortado_core.subprocess_discovery.concurrency_trees', 'cortado_core.subprocess_discovery.subtree_mining.right_most_path_extension',
              'cortado_core.subprocess_discovery.subtree_mining.metrics', 'cortado_core.subprocess_discovery.subtree_mining.blanket_mining', 'cortado_core.subprocess_discovery.subtree_mining.maximal_connected_components',
              'cortado_core.alignments.suffix_alignments'],
    install_requires=[
        'matplotlib==3.6.2',
        'cycler',
        'antlr4-python3-runtime==4.8',
        'networkx==3.0rc1',
        'pm4py==2.2.16',
        'ortools',
        'tqdm==4.64.1'
    ],
    python_requires=">=3.8"
)
