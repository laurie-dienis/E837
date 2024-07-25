**Analysis of experiment e837**

To convert the data : 

```
actroot -r tpc
```
```
actroot -r sil
```
```
actroot -f
```
```
actroot -m
```

To visualize the data event by event: 
```
actplot -v
```

To lauch the post-analysis code : 
```
cd PostAnalysis/
root 'Runner.cxx("number_of_the_pipeline")'
```
