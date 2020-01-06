RDKit Cookbook v2
%%%%%%%%%%%%%%%%%

.. contents:: :local:

Introduction
************

What is this?
=============

This document provides examples of how to carry out particular tasks using the RDKit functionality
from Python. The contents have been contributed by the RDKit community, tested with the latest 
RDKit release, and then compiled into this document. Note that this document is the second 
iteration of the Cookbook (i.e., v2). The old Cookbook written in Markdown is no longer 
maintained, but is available in prior RDKit releases for reference. The RDKit Cookbook v2 
is written in reStructuredText, which supports automatic Sphinx doctests, allowing for easier 
validation and maintenance of the RDKit Cookbook v2 code examples, where appropriate. 

What gets included?
===================

The examples included come from various online sources such as blogs, shared gists, and 
the RDKit mailing lists.  Generally, only minimal editing is added to the examples for 
formatting consistency. We have made a conscious effort to appropriately credit the original 
source and authors. For now, we are including anything that seems useful and reusable. 
One of the first priorities of this document is to compile useful examples shared on the RDKit 
mailing lists, as these can be difficult to discover. It will take some time, but we hope to expand 
this document into 100s of examples. As the document grows, it may make sense to prioritize 
examples included in the RDKit Cookbook v2 based on community demand. Another thought may be 
to turn the Cookbook v2 into a micro-publishing opportunity for the RDKit Community; that is, 
contributions to the RDKit Cookbook would be reviewed, and then assigned a unique DOI.

Feedback
========

If you have suggestions for how to improve the Cookbook v2 and/or examples you would like 
included, please contribute directly in the source document (the .rst file). The Index ID# 
is simply a way to track entries, new additions are sequentially numbered. Alternatively, 
you can also send Cookbook revision and and addition requests to the mailing list:
<rdkit-discuss@lists.sourceforge.net> (you will need to subscribe first).


Drawing Molecules (in a Jupyter Environment)
*******

Include an Atom Index
=====================

| **Author:** Takayuki Serizawa
| **Source:** `<https://iwatobipen.wordpress.com/2017/02/25/draw-molecule-with-atom-index-in-rdkit/>`_
| **Index ID#:** RDKitCB_0
| **Summary:** Draw a molecule with atom index numbers.

.. doctest::

   >>> from rdkit import Chem
   >>> from rdkit.Chem.Draw import IPythonConsole
   >>> from rdkit.Chem import Draw
   >>> IPythonConsole.ipython_useSVG=False
   >>> import rdkit
   >>> rdkit.__version__
   '2019.09.2'

.. doctest::
  
   >>> def mol_with_atom_index(mol):
   ...     atoms = mol.GetNumAtoms()
   ...     for idx in range(atoms):
   ...         mol.GetAtomWithIdx(idx).SetProp('molAtomMapNumber',str(mol.GetAtomWithIdx(idx).GetIdx()))
   ...     return mol

   >>> # Test in a kinase inhibitor
   >>> mol = Chem.MolFromSmiles("C1CC2=C3C(=CC=C2)C(=CN3C1)[C@H]4[C@@H](C(=O)NC4=O)C5=CNC6=CC=CC=C65")
   >>> # Default
   >>> mol # doctest: +ELLIPSIS
   <rdkit.Chem.rdchem.Mol object at 0x...>
   
.. image:: images/RDKitCB_0_im0.png

.. doctest::
  
   >>> # With atom index
   >>> mol_with_atom_index(mol) # doctest: +ELLIPSIS
   <rdkit.Chem.rdchem.Mol object at 0x...>

.. image:: images/RDKitCB_0_im1.png

Black and White Molecules
=====================

**Author:** Greg Landrum

**Source:** `<https://gist.github.com/greglandrum/d85d5693e57c306e30057ec4d4d11342>`_

**Index ID#:** RDKitCB_1

**Summary:** Draw a molecule in black and white.

.. doctest::

   >>> from rdkit import Chem
   >>> from rdkit.Chem.Draw import IPythonConsole
   >>> from rdkit.Chem import Draw
   >>> import rdkit
   >>> rdkit.__version__
   '2019.09.2'

.. doctest::

   >>> ms = [Chem.MolFromSmiles(x) for x in ('Cc1onc(-c2ccccc2)c1C(=O)N[C@@H]1C(=O)N2[C@@H](C(=O)O)C(C)(C)S[C@H]12','CC1(C)SC2C(NC(=O)Cc3ccccc3)C(=O)N2C1C(=O)O.[Na]')]
   >>> Draw.MolsToGridImage(ms) # doctest: +ELLIPSIS
   <PIL.PngImagePlugin.PngImageFile image mode=RGB size=600x200 at 0x...>

.. image:: images/RDKitCB_1_im0.png

.. doctest::

   >>> IPythonConsole.drawOptions.useBWAtomPalette()
   >>> Draw.MolsToGridImage(ms) # doctest: +ELLIPSIS
   <PIL.PngImagePlugin.PngImageFile image mode=RGB size=600x200 at 0x...>

.. image:: images/RDKitCB_1_im1.png

Highlight a Substructure in a Molecule
=====================

**Author:** Greg Landrum

**Source:** `<https://gist.github.com/greglandrum/5d45b56afe75603b955103cdd0d8e038>`_

**Index ID#:** RDKitCB_2

**Summary:** Draw a molecule with a substructure highlight.

.. doctest::

   >>> from rdkit import Chem
   >>> from rdkit.Chem.Draw import IPythonConsole
   >>> import rdkit
   >>> rdkit.__version__
   '2019.09.2'

.. doctest::

   >>> m = Chem.MolFromSmiles('c1cc(C(=O)O)c(OC(=O)C)cc1')
   >>> print(m.GetSubstructMatches(Chem.MolFromSmarts('C(=O)O')))
   ((3, 4, 5), (8, 9, 7))
   >>> m # doctest: +ELLIPSIS
   <rdkit.Chem.rdchem.Mol object at 0x...>

.. image:: images/RDKitCB_2_im0.png
   

Rings, Aromaticity, and Kekulization
************************************

Count Ring Systems
=====================

**Author:** Greg Landrum

**Source:** `<https://gist.github.com/greglandrum/de1751a42b3cae54011041dd67ae7415>`_

**Index ID#:** RDKitCB_3

**Summary:** Count ring systems in a molecule

.. doctest::

   >>> from rdkit import Chem
   >>> from rdkit.Chem.Draw import IPythonConsole

.. doctest::

   >>> def GetRingSystems(mol,includeSpiro=False):
   ...     ri = mol.GetRingInfo()
   ...     systems = []
   ...     for ring in ri.AtomRings():
   ...         ringAts = set(ring)
   ...         found = False
   ...         nSystems = []
   ...         for system in systems:
   ...             nInCommon = len(ringAts.intersection(system)) 
   ...             if nInCommon and (includeSpiro or nInCommon>1):
   ...                 ringAts = ringAts.union(system)
   ...             else:
   ...                 nSystems.append(system)
   ...         nSystems.append(ringAts)
   ...         systems = nSystems
   ...     return systems

.. doctest::

   >>> mol = Chem.MolFromSmiles('CN1C(=O)CN=C(C2=C1C=CC(=C2)Cl)C3=CC=CC=C3')
   >>> GetRingSystems(mol)
   [{1, 2, 4, 5, 6, 7, 8, 9, 10, 11, 12}, {14, 15, 16, 17, 18, 19}]

.. doctest::

   >>> # Draw molecule with atom index (see RDKitCB_0)
   >>> def mol_with_atom_index(mol):
   ...     atoms = mol.GetNumAtoms()
   ...     for idx in range(atoms):
   ...         mol.GetAtomWithIdx(idx).SetProp('molAtomMapNumber',str(mol.GetAtomWithIdx(idx).GetIdx()))
   ...     return mol
   >>> mol_with_atom_index(mol) # doctest: +ELLIPSIS
   <rdkit.Chem.rdchem.Mol object at 0x...>

.. image:: images/RDKitCB_3_im0.png

License
*******

.. image:: images/picture_5.png

This document is copyright (C) 2007-2020 by Greg Landrum and Vincent Scalfani.

This work is licensed under the Creative Commons Attribution-ShareAlike 4.0 License.
To view a copy of this license, visit http://creativecommons.org/licenses/by-sa/4.0/ 
or send a letter to Creative Commons, 543 Howard Street, 5th Floor, San Francisco, California, 94105, USA.


The intent of this license is similar to that of the RDKit itself. 
In simple words: “Do whatever you want with it, but please give us some credit.”
