#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright © 2019 yech <yech1990@gmail.com>
# Distributed under terms of the MIT license.
#
# Created: 2019-07-16 13:05

"""parse json format tree."""

import json

import pandas as pd


def parse_json_tree(
    json,
    name_label="name",
    children_label="children",
    length_label="branch_length",
):
    """ Return a newick string from a tree json format.
    Parameters:
    -----------
    json - A json object. See http://en.wikipedia.org/wiki/JSON
    Output:
    -------
    string - A string in newick format. See http://en.wikipedia.org/wiki/Newick_format
    Example:
    --------
    >>> json = {
        "name" : "F",
        "children": [
            {"name": "A", "branch_length": 0.1},
            {"name": "B", "branch_length": 0.2},
            {"name": "E","branch_length": 0.5,
            "children": [
                {"name": "C", "branch_length": 0.3},
                {"name": "D", "branch_length": 0.4}
                ]
            }
        ]
    }
    >>> print(json_to_newick(json))
    (A:0.1,B:0.2,(C:0.3,D:0.4)E:0.5)F;
    """

    record_list = []

    def _parse_json(json_obj):
        """Return a json object in the format desribed below
        (------- daughter node ------,------- daughter node ------, ...,------- daughter node ------)-- root/parent --
        ((children)name:branch_length,(children)name:branch_length, ...,(children)name:branch_length))name
        """

        record = {k: v for k, v in json_obj.items() if k != children_label}
        record["is_leaf"] = children_label not in json_obj
        record_list.append(record)
        record["short_name"] = record[name_label].split("_", 1)[-1]
        try:
            # Test is the key 'name' in the current level of the dictionary.
            newick = json_obj[name_label].split("_", 1)[-1]
        except KeyError:
            # Catch no 'name' trees and start the newick string with empty qoutes
            newick = ""

        if length_label in json_obj:
            newick = newick + ":" + str(json_obj[length_label])

        # If there are 'children'
        if children_label in json_obj:
            # Initialise a list to contain the daughter info
            info_list = []
            # For each child, treat it as a new dictionary object
            for child in json_obj[children_label]:
                # parse the newick string straight into it the list
                info_list.append(_parse_json(child))

            # join all the daughter info together with a comma
            info = ",".join(info_list)

            # Concatenate all the children together at the start of the parent newick string
            newick = "(" + info + ")" + newick

        return newick

    newick = _parse_json(json) + ";"
    table = pd.DataFrame(record_list)

    return newick, table


if __name__ == "__main__":
    with open(
        "./F6N3_sept2nd_2PART.withCells.allReads_cluster_colors.json", "r"
    ) as fp:
        obj = json.load(fp)[0]

    newick, table = parse_json_tree(obj)
    print(newick)
    table.to_csv("./test.csv", sep="\t", index=False)
