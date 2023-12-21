# Bushman Lab gene list

* this lab maintains a cancer gene list
* updated occasionally
* draws from several sources
* cited in several recent papers

```bash
wget -P data/ http://www.bushmanlab.org/assets/doc/allOnco_June2021.tsv

# Edit to add header entry for first column
sed -i '1s/^/"row"\t/' data/allOnco_June2021.tsv

./scripts/prepare.py 2>log.txt 1>prepared.tsv
```
