# Filogenomica_clase
Aquí explicare los pasos que realizamos en clase y analizaré los arboles que obtuvimos

## Pasos de lo que realizamos en clase
1. Primero entramos al cluster, a nuestras carpetas y creamos una carpeta para la clase y clonamos el repositorio con gitclone 
Gitclone link del repositorio y descomprimimos las carpetas con unzip
2.	abrimos salloc
3.	Activamos el ambiente con: ```conda activate biopt```
4.	Después de eso encontramos a los ortólogos con: ```orthofinder -os -M msa -S blast -f vertebrate_proteomes```
5.	Cambiamos el nombre de las especies con ```for f in *fa; do sed -E 's/_GENE.+//g' $f > out; mv out $f; done```, todo esto dentro del archivo ```Orthogroup_Sequences``` que se encontraba en la carpeta ```Orthogroups```
6.	Posteriormente hacemos un filtraje de calidad con ```for f in *fa; do prequal $f ; done```. Dentro del archivo ```Orthogroup_Sequences```. Todo esto con el fin de mejorar la calidad de nuestros resultados.
7.	Después realizamos el alineamiento multiple por medio de MAFFT ```for f in *filtered; do mafft $f > $f.mafft; done```.
En algunos casos este alineamiento puede ser erroneo porque sus secuencias son difíciles de alinear. A partir de esto removemos regiones muy incompletas. Para esto utilizamos BMGE y removimos posiciones que tengan gaps superiores al 80% 
```for f in *mafft; do bmge -i $f -t AA -g 0.8 -h 1 -w 1 -of $f.g08.fas; done```
8.	Posteriormente realizamos la concatenacion con FASconCAT 
```perl ./FASconCAT-G_v1.05.1.pl -l -s```
Ya teniendo mis secuencias concatenadas y que las consideramos como un solo locus que evoluciona de la misma forma, vamos a realizar dos árboles de máxima verosimilitud, uno más sofisticado que el otro.
9.	Para el primero, creamos un SBATCH con 10 núcleos y 20 minutos de tiempo con el siguiente comando
```module load iqtee```
```iqtree -s FcC_supermatrix.fas -m TEST -msub nuclear -bb 1000 -alrt 1000 -nt 1 -bnni -pre unpartitioned```
10.	Para el Segundo estamos delimitando donde inicia y donde termina cada gen para que calcule mejor el modelo para realizar las particiones y poder realizar el árbol. (PARA ESTE TAMBIÉN SE MANDA UN SBATCH CON LAS MISMAS CONDICIONES)

```module load iqtree```
```iqtree -s FcC_supermatrix.fas -spp FcC_supermatrix_partition.txt -m TEST -msub nuclear -merit AICc -bb 1000 -alrt 1000 -nt 1 -bnni -pre partitioned```

11.	Después de esto realizamos un análisis de coalescencia, donde utilizamos ASTRAL. Para ello creamos árboles de genes porque ASTRAL los toma como input.
```for f in *filtered.mafft.g08.fas; do iqtree -s $f -m TEST -msub nuclear -merit AICc -nt 1; done```
Posteriormente, colocamos los arboles en un mismo archivo
```cat *filtered.mafft.g08.fas.treefile > my_gene_trees.tre```
12.	Corrimos ASTRAL 
```java -jar /home/bio.pt/data/Astral/astral.5.7.8.jar -i my_gene_trees.tre -o species_tree_ASTRAL.tre 2> out.log```
13.	Para visualizar los arboles movimos los archivos al computador y los abrimos con FigTree.

## Resultados
# Partitioned tree
# Unpartitioned tree
# ASTRAL tree
