# Release process

Bump version

```bash
pip install bumpver
bumpver update

# For multiple releases within the same month
# bumpver update --patch
```

Build release assets

```bash
./scripts/build_release_assets.sh
```

After pushing new commits and tag, create release and upload build assets from `./build/<VERSION>/`
