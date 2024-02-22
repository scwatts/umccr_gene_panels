# Compiling Gene Utils from hmftools

Get hmftools GH repo and checkout target commit

```bash
git clone https://github.com/hartwigmedical/hmftools
cd hmftools/
git checkout a156ed6
```

Apply the following patch, `gene_utils.patch`

```patch
diff --git a/gene-utils/src/main/java/com/hartwig/hmftools/geneutils/drivers/GenerateDriverGeneFiles.java b/gene-utils/src/main/java/com/hartwig/hmftools/geneutils/drivers/GenerateDriverGeneFiles.java
index 647ea5525e..d9eebd6538 100644
--- a/gene-utils/src/main/java/com/hartwig/hmftools/geneutils/drivers/GenerateDriverGeneFiles.java
+++ b/gene-utils/src/main/java/com/hartwig/hmftools/geneutils/drivers/GenerateDriverGeneFiles.java
@@ -85,7 +85,6 @@ public class GenerateDriverGeneFiles
         List<DriverGene> driverGenes = DriverGeneFile.read(mDriverGenePanelFile);
         GU_LOGGER.info("loaded {} driver genes from {}", driverGenes.size(), mDriverGenePanelFile);

-        process(RefGenomeVersion.V37, driverGenes);
         process(RefGenomeVersion.V38, driverGenes);

         GU_LOGGER.info("file generation complete");
diff --git a/gene-utils/src/main/java/com/hartwig/hmftools/geneutils/fusion/GenerateFusionFiles.java b/gene-utils/src/main/java/com/hartwig/hmftools/geneutils/fusion/GenerateFusionFiles.java
index c88a1a3a5c..64bca80083 100644
--- a/gene-utils/src/main/java/com/hartwig/hmftools/geneutils/fusion/GenerateFusionFiles.java
+++ b/gene-utils/src/main/java/com/hartwig/hmftools/geneutils/fusion/GenerateFusionFiles.java
@@ -84,7 +84,6 @@ public class GenerateFusionFiles
             }
         }

-        createFusionFiles(RefGenomeVersion.V37, fusionRefData);
         createFusionFiles(RefGenomeVersion.V38, fusionRefData);

         GU_LOGGER.info("fusion reference file generation complete");
```

Applying...

```bash
patch -lp1 < gene_utils.patch
```

Build within Docker image

```bash
docker run -ti -v $(pwd):$(pwd) -w $(pwd) ubuntu:23.10

# NOTE(SW): MySQL and hmfpatients_test database are required to pass hmf-common tests
apt-get update && apt-get install -y maven mysql-server

service mysql start

mysql -u root -t <<EOF
CREATE DATABASE hmfpatients_test;
CREATE USER 'build'@'localhost' IDENTIFIED BY 'build';
GRANT ALL PRIVILEGES ON *.* TO 'build'@'localhost';
EOF

mvn install -pl gene-utils -am
```

Move JAR to `other/gene-utils_a156ed6.jar`
