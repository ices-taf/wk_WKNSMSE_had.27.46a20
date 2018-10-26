# wk_WKNSMSE_had.27.46a20

*note!*

This repo contains submodules.  To update the submodules to the *commit referenced in the github repository*  run:
```
git submodule update
```

And to *update the reference* to the most recent commit in the tools repo, go into that directory and perform a git pull
```
cd wk_WKNSMSE_tools
git pull origin master
```

and then navigate back to the repository root directory and commit and push the update
```
cd ..
git commit -am"updated submodule refs
git push
```


Some useful links:
* https://git-scm.com/book/en/v2/Git-Tools-Submodules
* https://blog.github.com/2016-02-01-working-with-submodules/
