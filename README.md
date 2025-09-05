
<b><p align='center'>[![Packt Sale](https://static.packt-cdn.com/assets/images/packt+events/Improve_UX.png)](https://packt.link/algotradingpython)</p></b> 




# Extending and Modifying LAMMPS Writing Your Own Source Code

<a href="https://www.packtpub.com/programming/extending-and-modifying-lammps-writing-your-own-source-code?utm_source=github&utm_medium=repository&utm_campaign=9781800562264"><img src="https://static.packt-cdn.com/products/9781800562264/cover/smaller" alt="Extending and Modifying LAMMPS Writing Your Own Source Code" height="256px" align="right"></a>

This is the code repository for [Extending and Modifying LAMMPS Writing Your Own Source Code](https://www.packtpub.com/programming/extending-and-modifying-lammps-writing-your-own-source-code?utm_source=github&utm_medium=repository&utm_campaign=9781800562264), published by Packt.

**A pragmatic guide to extending LAMMPS as per custom simulation requirements**

## What is this book about?
LAMMPS is one of the most widely used tools for running simulations for research in molecular dynamics. While the tool itself is fairly easy to use, more often than not you'll need to customize it to meet your specific simulation requirements. Extending and Modifying LAMMPS bridges this learning gap and helps you achieve this by writing custom code to add new features to LAMMPS source code. Written by ardent supporters of LAMMPS, this practical guide will enable you to extend the capabilities of LAMMPS with the help of step-by-step explanations of essential concepts, practical examples, and self-assessment questions. 

This book covers the following exciting features:
* Identify how LAMMPS input script commands are parsed within the source code
* Understand the architecture of the source code
* Relate source code elements to simulated quantities
* Learn how stored quantities are accessed within the source code
* Explore the mechanisms controlling pair styles, computes, and fixes
* Modify the source code to implement custom features in LAMMPS

If you feel this book is for you, get your [copy](https://www.amazon.com/dp/1800562268) today!

<a href="https://www.packtpub.com/?utm_source=github&utm_medium=banner&utm_campaign=GitHubBanner"><img src="https://raw.githubusercontent.com/PacktPublishing/GitHub/master/GitHub.png" 
alt="https://www.packtpub.com/" border="5" /></a>

## Instructions and Navigations
All of the code is organized into folders. For example, Chapter02.

The code will look like the following:
```
units metal
dimension 3
boundary p p p
atom_style atomic
atom_modify map array
```
The LAMMPS stable version 3Mar20 used in this book can be accessed from https://github.com/lammps/lammps/tree/stable_3Mar2020, or the zipped version can be downloaded from https://lammps.sandia.gov/tars/. The custom codes written in Chapters 9, 10, and 11, and Appendices B, C and D can be accessed from the various folders in this page labeled accordingly.

**Following is what you need for this book:**
This book is for students, faculty members, and researchers who are currently using LAMMPS or considering switching to LAMMPS, have a basic knowledge of how to use LAMMPS, and are looking to extend LAMMPS source code for research purposes. This book is not a tutorial on using LAMMPS or writing LAMMPS scripts, and it is assumed that the reader is comfortable with the basic LAMMPS syntax. The book is geared toward users with little to no experience in source code editing. Familiarity with C++ programming is helpful but not necessary.

With the following software and hardware list you can run all code files present in the book (Chapter 1-15).
### Software and Hardware List
| Chapter | Software required | OS required |
| -------- | ------------------------------------ | ----------------------------------- |
| 1-15 | LAMMPS stable version 3Mar20 | Windows, Mac OS X, and Linux (Any) |

We also provide a PDF file that has color images of the screenshots/diagrams used in this book. [Click here to download it]( https://static.packt-cdn.com/downloads/9781800562264_ColorImages.pdf).

### Related products
* Modern C++ Programming Cookbook - Second Edition [[Packt]](https://www.packtpub.com/product/modern-c-programming-cookbook-second-edition/9781800208988?utm_source=github&utm_medium=repository&utm_campaign=9781800208988) [[Amazon]](https://www.amazon.com/dp/1800208987)

* Windows Subsystem for Linux 2 (WSL 2) Tips, Tricks, and Techniques [[Packt]](https://www.packtpub.com/product/windows-subsystem-for-linux-2-wsl-2-tips-tricks-and-techniques/9781800562448?utm_source=github&utm_medium=repository&utm_campaign=9781800562448) [[Amazon]](https://www.amazon.com/dp/1800562446)

## Get to Know the Author
**Dr. Shafat Mubin**
(PhD, Physics, Penn State) is an assistant professor of physics at Valdosta State University. Since his graduate student days, he has worked with molecular simulations using primarily LAMMPS and has investigated a variety of simulation systems employing a wide array of techniques. He possesses extensive experience in writing custom routines and extending the LAMMPS source code, and hosts his own website to instruct and demonstrate the same to other users. At present, he is engaged in computational physics research including molecular simulations, and endeavours to train undergraduate students in computational techniques to help them better prepare for careers in physics.

**Jichen Li**
(graduated from Qingdao University of Science and Technology) is now studying for his master's degree at the University of Science and Technology of China. He used LAMMPS to conduct many molecular simulations to explore the relationship between polymer microstructure and macro mechanical and rheological properties. He developed several modeling and post-processing frameworks for LAMMPS and had a certain understanding of its program architecture. He dedicated to the community construction and wrote many columns and tutorials for LAMMPS starters. At present, he is working on the trans-scale simulation and the combination of deep learning and simulation.
### Download a free PDF

 <i>If you have already purchased a print or Kindle version of this book, you can get a DRM-free PDF version at no cost.<br>Simply click on the link to claim your free PDF.</i>
<p align="center"> <a href="https://packt.link/free-ebook/9781800562264">https://packt.link/free-ebook/9781800562264 </a> </p>