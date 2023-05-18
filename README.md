# Filogenomica_clase
Aquí explicare los pasos que realizamos en clase y analizaré los arboles que obtuvimos

## Pasos de lo que realizamos en clase
1. Primero entramos al cluster, a nuestras carpetas y creamos una carpeta para la clase y clonamos el repositorio con gitclone 
```Gitclone link del repositorio``` y descomprimimos las carpetas con unzip
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
### Partitioned tree
![image](https://github.com/MarianaT12/Filogenomicaclase/assets/130588298/c8db1027-503f-44ee-994e-7a202907cba6)
### Unpartitioned tree
![image](https://github.com/MarianaT12/Filogenomicaclase/assets/130588298/0e6ab5f7-225e-49e9-92f0-a2e4df5412e7)

* **Análisis**

Lo que se puede observar en ambos arboles son dos clados, donde el clado que se encuentra en la parte superior se considera monofilético y agrupa los peces óseos de la clase Actinopterigii. Ahora el clado de abajo que contiene especies de la clase Sarcopterigii se puede considerar parafilético porque se está excluyendo a _Latimeria chalumnae_, si lo hablamos desde el punto de vista sensu stricto (solo se tienen en cuenta las especies que se encuentran más estrechamente relacionadas). En el caso sensu lato, se consideraría monofilético (se tienen en cuenta todos los grupos que se derivan de un ancestro en común en los peces de aletas lobuladas), donde se observan los grupos de reptiles y mamíferos en dos grupos bien definidos.

En general, en ambos arboles las relaciones son las mismas y las diferencias en el numero de sustituciones son similares (longitud de las ramas). Sin embargo, los valores de los bootstraps son los que cambian en mayor medida, donde los valores del árbol realizado con las particiones es el que contiene mayor soporte. La diferencia principal radica en que en el modelo realizado con las particiones nosotros señalamos donde inicia y donde termina cada gen, reconociendo que los fragmentos de las secuencias moleculares pueden evolucionar de manera distinta, pues se considera que cada partición puede tener su propio modelo de sustitución, tasas de evolución y/o parámetros específicos. En contraste con el modelo que se realiza sin particiones donde a todas las secciones de la secuencia se tratan como una sola partición, sin tener en cuenta las diferencias potenciales en las tasas de evolución entre las diferentes de la secuencia.

### ASTRAL tree
![image](https://github.com/MarianaT12/Filogenomicaclase/assets/130588298/4e545e7d-0257-41bd-af55-8e9908a19a44)

* **Análisis**

Se observan 3 clados en el árbol filogenético (enumerados de arriba hacia debajo de forma ascendente), donde el clado de la parte superior es monofilético y se pueden observar el grupo de los reptiles y el de los mamíferos bien diferenciados. Ahora, en el caso del grupo de la mitad, sucede lo mismo que en los arboles de máxima verosimilitud. Si se observa desde el punto de vista sensu stricto, puede ser una relación parafilética porque en la clase Sarcopterygii se está excluyendo una especie (_Latimeria chalumnae_), sin embargo, desde el sensu lato se evidencia una relación monofilética. Ahora en el tercer clado, también se puede considerar monofilético desde el sensu lato, sin embargo si tenemos en cuenta las clases en sentido sensu stricto, podría ser parafiletico porque entre los Actinopterygii incluye un Sacopterygii y un Chondrichthyes como ancestros.

Observamos un soporte menor en el árbol de ASTRAL con respecto a los arboles de máxima verosimilitud, con algunas excepciones. ASTRAL sigue un enfoque coalescente, los cuales son capaces de capturar eventos de reticulación (interacción entre diferentes líneas de descendencia o linajes que resultan en un patrón de ramificación más complejo) o hibridación (proceso en el cual dos individuos se cruzan bien sean de especies estrechamente relacionadas o más distantes) en la historia evolutiva, mientras que los árboles de máxima verosimilitud asumen una evolución bifurcada sin eventos de intercambio genético entre las especies, ósea que si hay evidencia de alguno de estos dos eventos, es probable que el árbol muestre relaciones filogenéticas diferentes a las de los árboles de máxima verosimilitud. 


