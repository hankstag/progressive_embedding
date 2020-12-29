# Progressive Embedding

![](figure/teaser.png)

[Hanxiao Shen](http://cs.nyu.edu/~hanxiao/) [Zhongshi Jiang](http://cs.nyu.edu/~zhongshi/), [Denis Zorin](https://cims.nyu.edu/gcl/denis.html), [Daniele Panozzo](http://cs.nyu.edu/~panozzo/)<br/>
*ACM Transaction on Graphics (Proceedings of SIGGRAPH 2019)*<br/>
DOI: 10.1145/3306346.3323012

## Abstract
Tutte embedding is one of the most common building blocks in geome- try processing algorithms due to its simplicity and provable guarantees. Although provably correct in infinite precision arithmetic, it fails in chal- lenging cases when implemented using floating point arithmetic, largely due to the induced exponential area changes.
We propose Progressive Embedding, with similar theoretical guarantees to Tutte embedding, but more resilient to the rounding error of floating point arithmetic. Inspired by progressive meshes, we collapse edges on an invalid embedding to a valid, simplified mesh, then insert points back while maintaining validity. We demonstrate the robustness of our method by computing embeddings for a large collection of disk topology meshes.
By combining our robust embedding with a variant of the matchmaker algorithm, we propose a general algorithm for the problem of mapping multiply connected domains with arbitrary hard constraints to the plane, with applications in texture mapping and remeshing.

## To compile

Use the following cmds:
```
mkdir build
cd build
cmake -DCMAKE_BUILD_TYPE=Release ..
make
```

## To run

There are 4 executable files genearted as the compilation finishes:

- untangle_bin
- genus_zero_tutte_bin
- random_init_bin
- matchmaker_bin

---

Generate a failed case on Tutte's Embedding and then use it as initialization to generate a non-flip mapping with Progressive Embedding:
```
mkdir output
./genus_zero_tutte_bin --in ../data/62415_sf.obj -o ../data/62415_tutte_fail.obj
./untangle_bin --in ../data/62415_tutte_fail.obj -o output/62415_no_flip.obj
```

---

Genearte a random initialization with convex boundary and then use Progressive Embedding to fix the flips:
```
./random_init_bin --in ../data/retinal_miq.obj
./untangle_bin --in ../data/retinal_miq.obj_rand.obj -e 1
```

---

Randomly pick three vertices on mesh, pin them on the 2D plane, then generate a map
strictly satisfying the constraints using Matchmaker++:
```
./matchmaker_bin --in ../data/camel_miq.obj
```

## Data

To download the result data in our paper, please checkout this [link](https://drive.google.com/file/d/1caGIzPd9trlx0EvBbE06S2L3g1kHvIwJ/view?usp=sharing).

