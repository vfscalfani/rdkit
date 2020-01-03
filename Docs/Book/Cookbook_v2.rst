RDKit Cookbook v2
%%%%%%%%%%%%%%%%%

.. contents:: :local:

Introduction
************

What is this?
=============

This document provides examples of how to carry out particular tasks using the RDKit functionality from Python. The contents have been contributed by the RDKit community, tested with the latest RDKit release, and then compiled into this document. Note that this document is the second iteration of the Cookbook (i.e., v2). The old Cookbook written in markdown is no longer maintained, but is available in prior RDKit releases for reference. The RDKit Cookbook v2 is written in reStructuredText, which supports Sphinx doctests, allowing for easier validation and maintenance of the RDKit Cookbook v2 code examples, where appropiate. 

What gets included?
===================

The examples included come from various online sources such as blogs, shared gists, and the RDKit mailing lists. Generally, only minimal editing is added to the examples for formatting consistency. We have made a conscious effort to appropriately credit the original source and authors. For now, we are including anything that seems useful and reusable. Eventually, it may make sense to priortize examples included in the RDKit Cookbook v2 based on community demand. Another thought may be to turn the Cookbook v2 into a micropublishing oppurtnity for the RDKit Community; that is, contributions to the RDKit Cookbook would be reviewed, and then assigned a unique DOI. 

Feedback
========

If you have suggestions for how to improve the Cookbook and/or examples you would like included, please contribute directly in the source document (the .rst file) or send them to the mailing list: <rdkit-discuss@lists.sourceforge.net> (you will need to subscribe first).


Drawing Molecules (in a Jupyter Environment)
*******

Include an Atom Index
=========

**Author:** Takayuki Serizawa

**Source:** `<https://iwatobipen.wordpress.com/2017/02/25/draw-molecule-with-atom-index-in-rdkit/>`_

**Index ID#:** RDKitCB_0

**Summary:** Draw a molecule with index numbers.

   >>> from rdkit import Chem
   >>> from rdkit.Chem.Draw import IPythonConsole
   >>> from rdkit.Chem import Draw
   >>> IPythonConsole.ipython_useSVG=False
   >>> import rdkit
   >>> rdkit.__version__
   '2019.09.2'

   >>> def mol_with_atom_index(mol):
           atoms = mol.GetNumAtoms()
           for idx in range(atoms):
               mol.GetAtomWithIdx(idx).SetProp('molAtomMapNumber',str(mol.GetAtomWithIdx(idx).GetIdx()))
           return mol
   
   >>> mol = Chem.MolFromSmiles("C1CC2=C3C(=CC=C2)C(=CN3C1)[C@H]4[C@@H](C(=O)NC4=O)C5=CNC6=CC=CC=C65") # Test in a kinase inhibitor
   >>> mol # Default

.. image:: images/RDKitCB_0_im0.png

   >>> mol_with_atom_index(mol) # With index

.. image:: images/RDKitCB_0_im1.png

License
*******

.. image:: images/picture_5.png

This document is copyright (C) 2007-2020 by Greg Landrum and Vincent Scalfani.

This work is licensed under the Creative Commons Attribution-ShareAlike 4.0 License.
To view a copy of this license, visit http://creativecommons.org/licenses/by-sa/4.0/ or send a letter to Creative Commons, 543 Howard Street, 5th Floor, San Francisco, California, 94105, USA.


The intent of this license is similar to that of the RDKit itself. In simple words: “Do whatever you want with it, but please give us some credit.”
