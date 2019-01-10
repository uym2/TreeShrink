#!/bin/bash

vimdiff <(nl -ba $1) <(nl -ba $2)
