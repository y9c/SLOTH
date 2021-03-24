#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright © 2019 yech <yech1990@gmail.com>
# Distributed under terms of the MIT license.
#
# Created: 2019-05-22 20:27

"""shared function of this project."""

from __future__ import annotations

import logging

## global variable

_SEQ_COMPLEMENTS = str.maketrans("ACTG", "TGAC")


def get_logger(name: str) -> logging.Logger:
    """global logging."""
    logger: logging.Logger = logging.getLogger(name)
    if not logger.handlers:
        handler: logging.StreamHandler = logging.StreamHandler()
        formatter: logging.Formatter = logging.Formatter(
            "%(asctime)s %(name)-12s %(levelname)-8s %(message)s"
        )
        handler.setFormatter(formatter)
        logger.addHandler(handler)
        logger.setLevel(logging.DEBUG)
        #  logger.setLevel(logging.INFO)
    return logger


LOGGER: logging.Logger = get_logger(__name__)


def seq_revcom(seq):
    return seq.translate(_SEQ_COMPLEMENTS)[::-1]
