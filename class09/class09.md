class09
================

Preparing Data
--------------

``` r
url <- "https://bioboot.github.io/bimm143_S18/class-material/WisconsinCancer.csv"

wisc.df <- read.csv(url)

#get a matrix of column 3-32
wisc.data <- as.matrix(wisc.df[,3:32])

# Set the row names of wisc.data
row.names(wisc.data) <- wisc.df$id
#head(wisc.data) this views first couple
```

``` r
# Create diagnosis vector by completing the missing code
diagnosis <- as.numeric(wisc.df$diagnosis=="M")

#check how many should be
table(wisc.df$diagnosis == "M")
```

    ## 
    ## FALSE  TRUE 
    ##   357   212

``` r
#check diagnosis gives correct value
sum(diagnosis)
```

    ## [1] 212

Q1. How many observations are in this dataset?

``` r
nrow(wisc.df)
```

    ## [1] 569

Q2. How many variables/features in the data are suffixed with \_mean?

``` r
#show which values have _mean in its name 
grep("_mean", colnames(wisc.data))
```

    ##  [1]  1  2  3  4  5  6  7  8  9 10

``` r
#to count how many there are
length(grep("_mean", colnames(wisc.data)))
```

    ## [1] 10

Q3. How many of the observations have a malignant diagnosis?

``` r
table(wisc.df$diagnosis == "M")
```

    ## 
    ## FALSE  TRUE 
    ##   357   212

PCA
---

``` r
#might need to scale data so they are comparable
#check means of data
colMeans(wisc.data)
```

    ##             radius_mean            texture_mean          perimeter_mean 
    ##            1.412729e+01            1.928965e+01            9.196903e+01 
    ##               area_mean         smoothness_mean        compactness_mean 
    ##            6.548891e+02            9.636028e-02            1.043410e-01 
    ##          concavity_mean     concave.points_mean           symmetry_mean 
    ##            8.879932e-02            4.891915e-02            1.811619e-01 
    ##  fractal_dimension_mean               radius_se              texture_se 
    ##            6.279761e-02            4.051721e-01            1.216853e+00 
    ##            perimeter_se                 area_se           smoothness_se 
    ##            2.866059e+00            4.033708e+01            7.040979e-03 
    ##          compactness_se            concavity_se       concave.points_se 
    ##            2.547814e-02            3.189372e-02            1.179614e-02 
    ##             symmetry_se    fractal_dimension_se            radius_worst 
    ##            2.054230e-02            3.794904e-03            1.626919e+01 
    ##           texture_worst         perimeter_worst              area_worst 
    ##            2.567722e+01            1.072612e+02            8.805831e+02 
    ##        smoothness_worst       compactness_worst         concavity_worst 
    ##            1.323686e-01            2.542650e-01            2.721885e-01 
    ##    concave.points_worst          symmetry_worst fractal_dimension_worst 
    ##            1.146062e-01            2.900756e-01            8.394582e-02

``` r
#check standard deviation
apply(wisc.data,2,sd)
```

    ##             radius_mean            texture_mean          perimeter_mean 
    ##            3.524049e+00            4.301036e+00            2.429898e+01 
    ##               area_mean         smoothness_mean        compactness_mean 
    ##            3.519141e+02            1.406413e-02            5.281276e-02 
    ##          concavity_mean     concave.points_mean           symmetry_mean 
    ##            7.971981e-02            3.880284e-02            2.741428e-02 
    ##  fractal_dimension_mean               radius_se              texture_se 
    ##            7.060363e-03            2.773127e-01            5.516484e-01 
    ##            perimeter_se                 area_se           smoothness_se 
    ##            2.021855e+00            4.549101e+01            3.002518e-03 
    ##          compactness_se            concavity_se       concave.points_se 
    ##            1.790818e-02            3.018606e-02            6.170285e-03 
    ##             symmetry_se    fractal_dimension_se            radius_worst 
    ##            8.266372e-03            2.646071e-03            4.833242e+00 
    ##           texture_worst         perimeter_worst              area_worst 
    ##            6.146258e+00            3.360254e+01            5.693570e+02 
    ##        smoothness_worst       compactness_worst         concavity_worst 
    ##            2.283243e-02            1.573365e-01            2.086243e-01 
    ##    concave.points_worst          symmetry_worst fractal_dimension_worst 
    ##            6.573234e-02            6.186747e-02            1.806127e-02

``` r
# Perform PCA on wisc.data by completing the following code
wisc.pr <- prcomp( wisc.data, scale=T )
summary(wisc.pr)
```

    ## Importance of components:
    ##                           PC1    PC2     PC3     PC4     PC5     PC6
    ## Standard deviation     3.6444 2.3857 1.67867 1.40735 1.28403 1.09880
    ## Proportion of Variance 0.4427 0.1897 0.09393 0.06602 0.05496 0.04025
    ## Cumulative Proportion  0.4427 0.6324 0.72636 0.79239 0.84734 0.88759
    ##                            PC7     PC8    PC9    PC10   PC11    PC12
    ## Standard deviation     0.82172 0.69037 0.6457 0.59219 0.5421 0.51104
    ## Proportion of Variance 0.02251 0.01589 0.0139 0.01169 0.0098 0.00871
    ## Cumulative Proportion  0.91010 0.92598 0.9399 0.95157 0.9614 0.97007
    ##                           PC13    PC14    PC15    PC16    PC17    PC18
    ## Standard deviation     0.49128 0.39624 0.30681 0.28260 0.24372 0.22939
    ## Proportion of Variance 0.00805 0.00523 0.00314 0.00266 0.00198 0.00175
    ## Cumulative Proportion  0.97812 0.98335 0.98649 0.98915 0.99113 0.99288
    ##                           PC19    PC20   PC21    PC22    PC23   PC24
    ## Standard deviation     0.22244 0.17652 0.1731 0.16565 0.15602 0.1344
    ## Proportion of Variance 0.00165 0.00104 0.0010 0.00091 0.00081 0.0006
    ## Cumulative Proportion  0.99453 0.99557 0.9966 0.99749 0.99830 0.9989
    ##                           PC25    PC26    PC27    PC28    PC29    PC30
    ## Standard deviation     0.12442 0.09043 0.08307 0.03987 0.02736 0.01153
    ## Proportion of Variance 0.00052 0.00027 0.00023 0.00005 0.00002 0.00000
    ## Cumulative Proportion  0.99942 0.99969 0.99992 0.99997 1.00000 1.00000

Q4. From your results, what proportion of the original variance is captured by the first principal components (PC1)? Q5. How many principal components (PCs) are required to describe at least 70% of the original variance in the data? Q6. How many principal components (PCs) are required to describe at least 90% of the original variance in the data?

``` r
#find what to plot to show the cumulative proportion
attributes(wisc.pr)
```

    ## $names
    ## [1] "sdev"     "rotation" "center"   "scale"    "x"       
    ## 
    ## $class
    ## [1] "prcomp"

``` r
#this shows its "x"

#add +1 to col because F gives a Zero value that shows up white
plot( wisc.pr$x[,1], wisc.pr$x[,2], col = diagnosis+1, 
     xlab = "PC1", ylab = "PC2")
```

![](class09_files/figure-markdown_github/unnamed-chunk-8-1.png)

``` r
# Repeat for components 1 and 3
plot(wisc.pr$x[, c(1, 3)], col = (diagnosis + 1), 
     xlab = "PC1", ylab = "PC3")
```

![](class09_files/figure-markdown_github/unnamed-chunk-9-1.png)

``` r
#Calculate the variance of each principal component by squaring the sdev component
pr.var <- wisc.pr$sdev^2

# Variance explained by each principal component: pve
pve <- pr.var / sum(pr.var)

plot(pve, xlab = "Principal Component", 
     ylab = "Proportion of Variance Explained", 
     ylim = c(0, 1), type = "o")
```

![](class09_files/figure-markdown_github/unnamed-chunk-10-1.png)

``` r
# Alternative scree plot of the same data, note data driven y-axis
barplot(pve, ylab = "Precent of Variance Explained",
     names.arg=paste0("PC",1:length(pve)), las=2, axes = FALSE)
axis(2, at=pve, labels=round(pve,2)*100 )
```

![](class09_files/figure-markdown_github/unnamed-chunk-11-1.png)

``` r
# Plot cumulative proportion of variance explained
plot(cumsum(pve), xlab = "Principal Component", 
     ylab = "Cumulative Proportion of Variance Explained", 
     ylim = c(0, 1), type = "o")
```

![](class09_files/figure-markdown_github/unnamed-chunk-12-1.png)

``` r
#to get two graphs next to each other
par(mfrow=c(1,2))
plot(pve, xlab = "Principal Component", 
     ylab = "Proportion of Variance Explained", 
     ylim = c(0, 1), type = "o")
plot(cumsum(pve), xlab = "Principal Component", 
     ylab = "Cumulative Proportion of Variance Explained", 
     ylim = c(0, 1), type = "o")
```

![](class09_files/figure-markdown_github/unnamed-chunk-13-1.png)

Q9. For the first principal component, what is the component of the loading vector (i.e. wisc.pr$rotation\[,1\]) for the feature concave.points\_mean?

Q10. What is the minimum number of principal components required to explain 80% of the variance of the data? 4

Hierarchical clustering of case data
------------------------------------

``` r
# Scale the wisc.data data: data.scaled
data.scaled <- scale(wisc.data)

data.dist <- dist(data.scaled)

#Create a hierarchical clustering model using complete linkage
wisc.hclust <- hclust(data.dist, method = "complete")

hcluster <- plot(wisc.hclust)
```

![](class09_files/figure-markdown_github/unnamed-chunk-14-1.png)

to analyze this, cut it into different clusters

``` r
wisc.hclust.clusters <- cutree(wisc.hclust, k=4)

wisc.hclust.clusters
```

    ##    842302    842517  84300903  84348301  84358402    843786    844359 
    ##         1         1         1         2         1         1         1 
    ##  84458202    844981  84501001    845636  84610002    846226    846381 
    ##         1         1         2         3         1         1         3 
    ##  84667401  84799002    848406  84862001    849014   8510426   8510653 
    ##         1         1         3         1         1         3         3 
    ##   8510824   8511133    851509    852552    852631    852763    852781 
    ##         3         1         1         1         1         1         1 
    ##    852973    853201    853401    853612  85382601    854002    854039 
    ##         1         3         1         1         1         1         1 
    ##    854253    854268    854941    855133    855138    855167    855563 
    ##         1         1         3         3         1         3         1 
    ##    855625    856106  85638502    857010  85713702     85715    857155 
    ##         1         1         1         1         3         1         3 
    ##    857156    857343    857373    857374    857392    857438  85759902 
    ##         3         3         3         3         1         3         3 
    ##    857637    857793    857810    858477    858970    858981    858986 
    ##         1         1         3         3         3         3         1 
    ##    859196  85922302    859283    859464    859465    859471    859487 
    ##         3         1         1         3         3         2         3 
    ##    859575    859711    859717    859983   8610175   8610404   8610629 
    ##         1         3         1         1         3         3         3 
    ##   8610637   8610862   8610908    861103   8611161   8611555   8611792 
    ##         1         2         3         3         1         1         1 
    ##   8612080   8612399  86135501  86135502    861597    861598    861648 
    ##         3         1         3         1         3         1         3 
    ##    861799    861853    862009    862028     86208     86211    862261 
    ##         3         3         3         1         1         3         3 
    ##    862485    862548    862717    862722    862965    862980    862989 
    ##         3         3         3         3         3         3         3 
    ##    863030    863031    863270     86355    864018    864033     86408 
    ##         1         1         3         1         3         3         3 
    ##     86409    864292    864496    864685    864726    864729    864877 
    ##         3         3         3         3         3         1         1 
    ##    865128    865137     86517    865423    865432    865468     86561 
    ##         3         3         1         2         3         3         3 
    ##    866083    866203    866458    866674    866714      8670  86730502 
    ##         1         3         1         1         3         1         1 
    ##    867387    867739    868202    868223    868682    868826    868871 
    ##         3         1         3         3         3         1         3 
    ##    868999    869104    869218    869224    869254    869476    869691 
    ##         3         3         3         3         3         3         1 
    ##  86973701  86973702    869931 871001501 871001502   8710441     87106 
    ##         3         3         3         3         3         2         3 
    ##   8711002   8711003   8711202   8711216    871122    871149   8711561 
    ##         3         3         1         3         3         3         3 
    ##   8711803    871201   8712064   8712289   8712291     87127   8712729 
    ##         1         1         3         1         3         3         3 
    ##   8712766   8712853  87139402     87163     87164    871641    871642 
    ##         1         3         3         3         1         3         3 
    ##    872113    872608  87281702    873357    873586    873592    873593 
    ##         3         3         1         3         3         1         1 
    ##    873701    873843    873885    874158    874217    874373    874662 
    ##         1         3         1         3         3         3         3 
    ##    874839    874858    875093    875099    875263  87556202    875878 
    ##         3         2         3         3         1         1         3 
    ##    875938    877159    877486    877500    877501    877989    878796 
    ##         1         3         1         1         3         3         1 
    ##     87880     87930    879523    879804    879830   8810158   8810436 
    ##         1         3         3         3         3         1         3 
    ## 881046502   8810528   8810703 881094802   8810955   8810987   8811523 
    ##         1         3         4         3         1         1         3 
    ##   8811779   8811842  88119002   8812816   8812818   8812844   8812877 
    ##         3         1         1         3         3         3         1 
    ##   8813129  88143502  88147101  88147102  88147202    881861    881972 
    ##         3         3         3         3         3         1         1 
    ##  88199202  88203002  88206102    882488  88249602  88299702    883263 
    ##         3         3         1         3         3         1         1 
    ##    883270  88330202  88350402    883539    883852  88411702    884180 
    ##         3         1         3         3         3         3         1 
    ##    884437    884448    884626  88466802    884689    884948  88518501 
    ##         3         3         3         3         3         1         3 
    ##    885429   8860702    886226    886452  88649001    886776    887181 
    ##         1         3         1         3         1         1         1 
    ##  88725602    887549    888264    888570    889403    889719  88995002 
    ##         1         1         3         1         3         1         1 
    ##   8910251   8910499   8910506   8910720   8910721   8910748   8910988 
    ##         3         3         3         3         3         3         1 
    ##   8910996   8911163   8911164   8911230   8911670   8911800   8911834 
    ##         3         3         3         3         3         3         3 
    ##   8912049   8912055     89122   8912280   8912284   8912521   8912909 
    ##         1         3         1         1         3         3         3 
    ##      8913   8913049  89143601  89143602      8915    891670    891703 
    ##         3         3         3         3         3         3         3 
    ##    891716    891923    891936    892189    892214    892399    892438 
    ##         3         3         3         3         3         3         1 
    ##    892604  89263202    892657     89296    893061     89344     89346 
    ##         3         1         3         3         3         3         3 
    ##    893526    893548    893783  89382601  89382602    893988    894047 
    ##         3         3         3         3         3         3         3 
    ##    894089    894090    894326    894329    894335    894604    894618 
    ##         3         3         1         3         3         3         3 
    ##    894855    895100  89511501  89511502     89524    895299   8953902 
    ##         3         1         3         3         3         3         1 
    ##    895633    896839    896864    897132    897137    897374  89742801 
    ##         1         1         1         3         3         3         1 
    ##    897604    897630    897880     89812     89813    898143     89827 
    ##         3         1         3         1         3         3         3 
    ##    898431  89864002    898677    898678     89869    898690    899147 
    ##         1         3         3         3         3         3         3 
    ##    899187    899667    899987   9010018    901011   9010258   9010259 
    ##         3         1         1         1         3         3         3 
    ##    901028   9010333 901034301 901034302    901041   9010598   9010872 
    ##         3         3         3         3         3         3         3 
    ##   9010877    901088   9011494   9011495   9011971   9012000   9012315 
    ##         3         1         1         3         1         1         1 
    ##   9012568   9012795    901288   9013005    901303    901315   9013579 
    ##         3         1         1         3         3         3         3 
    ##   9013594   9013838    901549    901836     90250     90251    902727 
    ##         3         1         3         3         3         3         3 
    ##     90291    902975    902976    903011     90312  90317302    903483 
    ##         3         3         3         3         1         3         3 
    ##    903507    903516    903554    903811  90401601  90401602    904302 
    ##         1         1         3         3         3         3         3 
    ##    904357  90439701    904647    904689      9047    904969    904971 
    ##         3         1         3         3         3         3         3 
    ##    905189    905190  90524101    905501    905502    905520    905539 
    ##         3         3         1         3         3         3         3 
    ##    905557    905680    905686    905978  90602302    906024    906290 
    ##         3         3         3         3         1         3         3 
    ##    906539    906564    906616    906878    907145    907367    907409 
    ##         3         1         3         3         3         3         3 
    ##     90745  90769601  90769602    907914    907915    908194    908445 
    ##         3         3         3         1         3         1         1 
    ##    908469    908489    908916    909220    909231    909410    909411 
    ##         3         1         3         3         3         3         3 
    ##    909445  90944601    909777   9110127   9110720   9110732   9110944 
    ##         3         3         3         3         3         1         3 
    ##    911150 911157302   9111596   9111805   9111843    911201    911202 
    ##         3         1         3         1         3         3         3 
    ##   9112085   9112366   9112367   9112594   9112712 911296201 911296202 
    ##         3         3         3         3         3         1         4 
    ##   9113156 911320501 911320502   9113239   9113455   9113514   9113538 
    ##         3         3         3         3         3         3         1 
    ##    911366   9113778   9113816    911384   9113846    911391    911408 
    ##         1         3         3         3         3         3         3 
    ##    911654    911673    911685    911916    912193     91227    912519 
    ##         3         3         3         1         3         3         3 
    ##    912558    912600    913063    913102    913505    913512    913535 
    ##         3         3         3         3         1         3         3 
    ##  91376701  91376702    914062    914101    914102    914333    914366 
    ##         3         3         1         3         3         3         1 
    ##    914580    914769     91485    914862     91504     91505    915143 
    ##         3         1         1         3         1         3         1 
    ##    915186    915276  91544001  91544002    915452    915460     91550 
    ##         3         3         3         3         3         1         3 
    ##    915664    915691    915940  91594602    916221    916799    916838 
    ##         3         1         3         3         3         1         1 
    ##    917062    917080    917092  91762702     91789    917896    917897 
    ##         3         3         3         1         3         3         3 
    ##     91805  91813701  91813702    918192    918465     91858  91903901 
    ##         3         1         3         3         3         3         3 
    ##  91903902  91930402    919537    919555  91979701    919812    921092 
    ##         3         1         3         1         3         1         3 
    ##    921362    921385    921386    921644    922296    922297    922576 
    ##         3         3         1         3         3         3         3 
    ##    922577    922840    923169    923465    923748    923780    924084 
    ##         3         3         3         3         3         3         3 
    ##    924342    924632    924934    924964    925236    925277    925291 
    ##         3         3         3         3         3         3         3 
    ##    925292    925311    925622    926125    926424    926682    926954 
    ##         3         3         1         1         1         1         3 
    ##    927241     92751 
    ##         1         3

``` r
table(wisc.hclust.clusters, diagnosis)
```

    ##                     diagnosis
    ## wisc.hclust.clusters   0   1
    ##                    1  12 165
    ##                    2   2   5
    ##                    3 343  40
    ##                    4   0   2

try again with 6 clusters

``` r
wisc.hclust.clusters6 <- cutree(wisc.hclust, k=6)

table(wisc.hclust.clusters6, diagnosis)
```

    ##                      diagnosis
    ## wisc.hclust.clusters6   0   1
    ##                     1  12 165
    ##                     2   0   5
    ##                     3 331  39
    ##                     4   2   0
    ##                     5  12   1
    ##                     6   0   2

with 3

``` r
wisc.hclust.clusters3 <- cutree(wisc.hclust, k=3)

table(wisc.hclust.clusters3, diagnosis)
```

    ##                      diagnosis
    ## wisc.hclust.clusters3   0   1
    ##                     1 355 205
    ##                     2   2   5
    ##                     3   0   2

Clustering on PCA results
-------------------------

``` r
## Use the distance along the first 7 PCs for clustering i.e. wisc.pr$x[, 1:7]
# need to input distance vector
d.pr <- dist(wisc.pr$x[, 1:7])
wisc.pr.hclust <- hclust(d.pr, method="complete")
plot(wisc.pr.hclust)
```

![](class09_files/figure-markdown_github/unnamed-chunk-18-1.png) check how well it did

``` r
wisc.pr.hclust.clusters <- cutree(wisc.pr.hclust, k=4)
table(wisc.pr.hclust.clusters, diagnosis)
```

    ##                        diagnosis
    ## wisc.pr.hclust.clusters   0   1
    ##                       1   5 113
    ##                       2 350  97
    ##                       3   2   0
    ##                       4   0   2

this does worse at seperating the diagnoses into clusters than the hieracrchical clusters

``` r
#to predict whether these 2 people would have cancer or not

url2 <- "https://tinyurl.com/new-samples-CSV"
new <- read.csv(url2)
npc <- predict(wisc.pr, newdata=new)
plot(wisc.pr$x[,1:2], col=diagnosis+1)
points(npc[,1], npc[,2], col="blue", pch=16)
```

![](class09_files/figure-markdown_github/unnamed-chunk-20-1.png)
