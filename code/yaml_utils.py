#!/usr/bin/env python

'''
Author: Prag Batra prag@stanford.edu

Purpose:
    Helper methods for parsing YAML configuration files (e.g. for STMP).
    

Explanation:

    

Example:

    

'''
import yaml
import yaml_keys
import re
import os
from collections import OrderedDict

"""Make PyYAML output an OrderedDict.
It will do so fine if you use yaml.dump(), but that generates ugly, 
non-standard YAML code.
To use yaml.safe_dump(), you need the following.
(Credit goes to http://blog.elsdoerfer.name/2012/07/26/make-pyyaml-output-an-ordereddict/ and https://gist.github.com/miracle2k/3184458#file-odict-py)
"""

def represent_odict(dump, tag, mapping, flow_style=None):
    """Like BaseRepresenter.represent_mapping, but does not issue the sort().
    """
    value = []
    node = yaml.MappingNode(tag, value, flow_style=flow_style)
    if dump.alias_key is not None:
        dump.represented_objects[dump.alias_key] = node
    best_style = True
    if hasattr(mapping, 'items'):
        mapping = mapping.items()
    for item_key, item_value in mapping:
        node_key = dump.represent_data(item_key)
        node_value = dump.represent_data(item_value)
        if not (isinstance(node_key, yaml.ScalarNode) and not node_key.style):
            best_style = False
        if not (isinstance(node_value, yaml.ScalarNode) and not node_value.style):
            best_style = False
        value.append((node_key, node_value))
    if flow_style is None:
        if dump.default_flow_style is not None:
            node.flow_style = dump.default_flow_style
        else:
            node.flow_style = best_style
    return node

yaml.SafeDumper.add_representer(OrderedDict,
    lambda dumper, value: represent_odict(dumper, u'tag:yaml.org,2002:map', value))


# works for both list and dictionary of columns
def convertColumns(cols, yaml_commands):
    actualCols = []
    yaml_datasets = get_datasets(yaml_commands)
    for col in cols:
        if(isinstance(col, list)):
            convertedColList = convertColumns(col, yaml_commands) # recursion
            actualCols.append(convertedColList)
            continue
        if('.' not in col): # nothing to do here
            actualCols.append(col)
            continue
        #else
        colComponents = col.split('.')
        dataset_yaml = yaml_datasets[colComponents[0]]
        dataset_annotation_name = dataset_yaml[yaml_keys.kDAnnotation]
        if(dataset_yaml[yaml_keys.kDCategory] == yaml_keys.kDCategoryTypeRegion):
            dataset_annotation_name += '_r'
        actualCols.append(dataset_annotation_name+'_'+colComponents[1])

    return actualCols


# helper function to load yaml datasets as an ordered dictionary as opposed to a regular (unordered) dictionary
def ordered_load(stream, Loader=yaml.Loader, object_pairs_hook=OrderedDict):
    class OrderedLoader(Loader):
        pass
    def construct_mapping(loader, node):
        loader.flatten_mapping(node)
        return object_pairs_hook(loader.construct_pairs(node))
    OrderedLoader.add_constructor(
        yaml.resolver.BaseResolver.DEFAULT_MAPPING_TAG,
        construct_mapping)
    return yaml.load(stream, OrderedLoader)

# loads YAML file and returns content as nested python dictionaries/arrays
# NOTE: generally, you should use parse_yaml_input_files below instead
def parse_yaml(loc, load_ordered = False):
    # Parse a YAML file which will instruct the annotation steps to be completed.
    with open(loc, "r") as stream:
        if(load_ordered):
            yaml_commands = ordered_load(stream, yaml.SafeLoader)
        else:
            yaml_commands = yaml.safe_load(stream)
        return yaml_commands

# combines yaml input files, parses, and outputs single python structure containing commands (nested dictionaries, etc.)
def parse_yaml_input_files(dataset_file, modules_file):
    config_cmds = {}
    dataset_cmds = parse_yaml(dataset_file, load_ordered=True) # info about each dataset (including defaults)
    dataset_default_cmds = dataset_cmds[yaml_keys.kDDefaults]
    del dataset_cmds[yaml_keys.kDDefaults]
    module_cmds = parse_yaml(modules_file)
    #combine info to yield final config commands
    config_cmds[yaml_keys.kModules] = module_cmds
    config_cmds[yaml_keys.kDatasets] = dataset_cmds # = datasets since we deleted defaults from here
    config_cmds[yaml_keys.kDatasetDefaults] = dataset_default_cmds
    return config_cmds

# exports existing yaml commands to output YAML files (modules.yml and datasets.yml), which can be used to rerun STMP with the same configuration.
def write_output_yaml_files(yaml_commands, output_dir):
    #write datasets.yml
    dataset_defaults = {}
    dataset_defaults[yaml_keys.kDDefaults] = yaml_commands[yaml_keys.kDatasetDefaults]
    datasets = get_datasets(yaml_commands)
    #create 1 big ordereddict containing defaults followed by per-dataset info
    datasets_and_defaults = OrderedDict(list(dataset_defaults.items()) + list(datasets.items()))
    datasets_out = open(os.path.join(output_dir, 'datasets.yml'), 'w')
    yaml.safe_dump(datasets_and_defaults, datasets_out, default_flow_style=False)
    datasets_out.close()
    #write modules.yml
    modules = yaml_commands[yaml_keys.kModules]
    modules_out = open(os.path.join(output_dir, 'modules.yml'), 'w')
    yaml.safe_dump(modules, modules_out, default_flow_style=False)
    modules_out.close()


# gets absolute path wrt location of this file (which should be in the same directory as stmp.py)
def get_abs_path(yaml_path):
    script_dir = os.path.dirname(os.path.realpath(__file__))
    if(not yaml_path.startswith('/')): # not an absolute path
        yaml_path = os.path.join(script_dir, yaml_path)
    return yaml_path

# splits a given set of threshold values for a column using the separator specified in the YAML tiering config
def split_multiple_col_thresholds(col_threshold_str, yaml_commands):
    col_threshold_separator = yaml_commands[yaml_keys.kModules][yaml_keys.kTiering][yaml_keys.kTColMultipleThresholdSeparator]
    return col_threshold_str.split(col_threshold_separator)

# gets ordered list of all datasets (as an OrderedDict)
def get_datasets(yaml_commands):
    return yaml_commands[yaml_keys.kDatasets]

def get_dataset_defaults(yaml_commands):
    return yaml_commands[yaml_keys.kDatasetDefaults]


def is_region_dataset(dataset_yaml_name, yaml_commands):
    return get_datasets(yaml_commands)[dataset_yaml_name][yaml_keys.kDCategory] == yaml_keys.kDCategoryTypeRegion

# converts from annotated dataset name to top-level dataset name in yaml file
def annotated_to_yaml_dataset_name(annotated_dataset_name, yaml_commands):
    datasets_yaml = get_datasets(yaml_commands)
    for dataset in datasets_yaml:
        if(annotated_dataset_name == datasets_yaml[dataset][yaml_keys.kDAnnotation] or (datasets_yaml[dataset][yaml_keys.kDCategory] == yaml_keys.kDCategoryTypeRegion and re.sub('_r$', '', annotated_dataset_name) == datasets_yaml[dataset][yaml_keys.kDAnnotation])):
            return dataset
    #else
    raise ValueError('Could not find dataset in YAML: ' + str(annotated_dataset_name))

