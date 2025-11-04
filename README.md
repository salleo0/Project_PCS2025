# **Programming and Scientific Computing Project: Geodetic Polyhedra**

This project involves the C++ implementation of a program to generate and analyze **Geodesic Polyhedra** and their duals, **Goldberg Polyhedra**, based on the project requirements for the *Programmazione e Calcolo Scientifico* course.

## **Project Structure**

The assignment is divided into two main parts:

* **Part I:** Implementation of **Class I** polyhedra and their duals.  
* **Part II:** Extension of the algorithm to support **Class II** polyhedra.

## **Part I: Core Requirements**

This part focuses on **Class I** geodesic polyhedra, which are defined by an input quadruple (p, q, b, c) where either $b=0, c \\ge 1$ or $c=0, b \\ge 1$.

### **1\. Polyhedron Data Structure**

You must define a data structure to represent all properties of a polyhedron, organized into "cells":

* **0D Cells (Vertices):** Each vertex must have a unique ID (starting from 0\) and its 3D coordinates (x, y, z).  
* **1D Cells (Edges):** Each edge must have a unique ID, and the IDs of its origin and end vertices.  
* **2D Cells (Faces):** Each face must have a unique ID, the count of its vertices and edges, and ordered lists of the IDs for its constituent vertices and edges.  
* **3D Cells (Polyhedron):** The polyhedron itself must have a unique ID, the total count of its vertices, edges, and faces, and lists of all their IDs.

### **2\. Polyhedra Generation**

The program must handle the following generation tasks based on the input (p, q, b, c):

* **Geodesic Polyhedra (Class I):** If the input p=3, the program must generate the corresponding Class I geodesic polyhedron.  
* **Goldberg Polyhedra (Class I):** If the input q=3, the program must generate the corresponding Class I Goldberg polyhedron (the dual of the geodesic solid).

### **3\. Sphere Projection**

All vertices of the generated polyhedra must be projected so that they lie on the surface of a unit sphere (radius 1\) centered at the origin.

### **4\. Output and Visualization**

* **Text Files:** The program must output 4 .txt files detailing the generated structure: Cell0Ds.txt, Cell1Ds.txt, Cell2Ds.txt, and Cell3Ds.txt.  
* **ParaView:** The program must be able to print the vertices and edges to be visualized in ParaView. (Visualization of faces or the full 3D cell is not required).

### **5\. Shortest Path Algorithm**

The program must also handle a 6-tuple input: (p, q, b, c, id\_vertex\_1, id\_vertex\_2).

* **Task:** Find the shortest path between the two specified vertex IDs. The graph is defined by the polyhedron's vertices (nodes) and its edges (1D cells).  
* **Screen Output:** Print the number of edges that make up the shortest path and the sum of their lengths.  
* **ParaView Output:** Highlight the path in ParaView. Vertices and edges on the path should have a property ShortPath \= 1, while all others have ShortPath \= 0\.

### **6\. Validation and Testing**

* **Input Validation:** You must always verify the correctness of the user's input.  
* **Unit Testing:** Use **Google Test** to verify the correct functioning of every logical unit defined in your code.

## **Part II: Extension Requirements**

This part involves modifying the algorithm from Part I to add support for **Class II** geodesic solids.

* **Task:** Enable the construction of solids where $b=c$ (and $b \\ge 1$).

## **Deliverables**

To receive a positive evaluation, the following must be submitted:

1. The complete **C++ source code** (hosted on this shared GitHub repository).  
2. The corresponding **UML documentation** that describes the logical units defined in the source code.  
3. A **presentation** in which your results are shown.