#!/usr/bin/env python

"""
TODO: Modify module doc.
"""

from __future__ import division

__author__ = "Shyue Ping Ong"
__copyright__ = "Copyright 2012, The Materials Project"
__version__ = "0.1"
__maintainer__ = "Shyue Ping Ong"
__email__ = "shyuep@gmail.com"
__date__ = "7/22/15"


from flask import render_template, request, redirect, url_for, make_response,\
    session, jsonify, Markup


from pymatgen.gui import app


@app.route('/', methods=['GET'])
def index():
    """Initialize a session for the current user, and render index.html."""
    # Create a state token to prevent request forgery.
    # Store it in the session for later validation.
    response = make_response(render_template('index.html'))
    response.headers['Content-Type'] = 'text/html'
    return response
