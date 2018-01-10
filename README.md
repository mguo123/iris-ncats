# iris-ncats: A Conversational Agent for Biological Data Science 

This repository contains the prototype for the Stanford Effort for the Biomedical Data Translator. The system we developed is able to answer the following questions:

* Given a condition, are any genetic conditions potentially protective and why?
* Given a drug and a condition, construct a pathway description that explains how the drug effects its action. Additionally, generate hypotheses of other therapeutic uses for the given drug.

These queries take in relevant arguments, use a combination of semantic similarity comparison methods, genetic network analysis, and/or enrichment calculations to generate a plausible hypotheses. These hypotheses are meant to act as a launch point for further research efforts.

This prototype is based on the IRIS, which is an open source conversational agent, originally to facilitate tasks in data science for non-programmers.  We tailored it to a biomedical context, allowing it to recognize key entities (drugs, genes, diseases), pull from various biomedical data sources, construct queries based on free text inputted by the user, and perform various analyses to produce the final answer or dossier.

The original IRIS system was created by the Bernstein lab at Stanford Computer Science, and the [iris-agents](https://github.com/Ejhfast/iris-agent) repository can be found here.

![interface](/interface_ncats.png)


## Getting Starting 

### Note: This is currently an alpha prototype release. Additions and bux fixes are on going.


In terms of software structure, this effort can be thought of in a couple forms:
1. as a standalone desktop Electron app (for OSX and Linux with future compatibility with Windows)
2. as a web application that is ported through a basic http-server
3. as a command-line interface for individuals who want more programmatic control (repository can be found [here](https://github.com/alavertu/ncats_altman))

The backend of Iris is written in Python while the front end is written in Javascript/CSS/HTML.  These are instructions to install and run iris-ncats in debugging mode that can be used for developers. We recommend using OSX for running the debugging mode, as the Windows and Linux versions are still in pre-release. A self-contained Electron app for OSX will be available shortly.

### Prerequisites for Developing

* [anaconda 3](https://conda.io/docs/install/quick.html). (Make sure to run `source ~/.bash_profile` after you have installed Anaconda, if it is not appearing in your path.)
* [python 2.7](https://www.python.org/download/releases/2.7/). This is in order to install the javascript components
* [xcode](https://developer.apple.com/xcode/downloads/)
* node.js: 
    * run `ruby -e "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/install)"`
    * `brew install node` (if you hit errors: go [here](http://blog.teamtreehouse.com/install-node-js-npm-mac))
* graphvix - to install you need homobrew, then run the following lines 
```
brew install graphviz
pip install pygraphviz --install-option="--include-path=/usr/include/graphviz" --install-option="--library-path=/usr/lib/graphviz/"
```


### Install and run the Python components:
```
    cd [PATH_TO_IRIS]/backend
    # create new conda env called iris with necessary packages
    conda create --name iris --file packages.txt
    # enter the conda env
    source activate iris
    # install the remaining pip packages
    pip install -r requirements.txt
```

### Install and run the Javascript components:
```
    cd [PATH_TO_IRIS]
    # If you do not already have webpack
    npm install webpack -g
    # Now installing iris javascript components
    npm install --python=python2.7
    # build JS app with webpack (can also run webpack --watch in seperate command window)
    webpack
```



## Using Iris (desktop dev mode)
In order to run Iris in development mode. 

1. Open two terminals, one that runs frontend launcher commands, the other your Python-based backend.
2. In the python terminal type:
```
    source activate iris 


    # run the backend application
    cd [PATH_TO_IRIS]/backend/app # enter backend path

    python app.py
```
3. Wait until your backend has finished booting up (should take between 30-40 seconds) you should see in your terminal
```
======== Running on http://0.0.0.0:8000 ========
(Press CTRL+C to quit)
```

4. In the other window, type:
```

    cd [PATH_TO_IRIS]

    # start electron (this will open the application automatically)
    npm start
```

## Using Iris on the web
Iris can run within a desktop electron app (currently in dev mode only, but will become a standalone app) and as a web application.

To run as a web application follow steps 1-3 for the local version to start up the Python backend. Then for the web version enter in:
```
    cd [PATH_TO_IRIS]

    http-server . [PORT_NUM] # default is 8080


```

Open a webbrowser to the indicated url (i.e. `localhost:PORT_NUM`) and Iris should be displayed. 

Notes if you are developping on AWS or a remote server:
* We recommend creating 2 screens to run the front and backend simultaneously. This can be done as follows
```
    # in your AWS instance
    cd iris_ncats # PATH TO REPOSITORY

    # Start a new screen session
    screen -S python
    # in screen
    source activate iris
    cd backend/app
    python app.py

    # Once, backed is loaded, enter Ctrl+A then Ctrl+D to exit screen

    # Start the html screen session
    screen -S html
    # in screen 
    http-server .


```
* Once your program is up and running, you need to forward both the port number entered for you http-server as well as for the python backend. This can be done via the ssh command.
```
    # on your local computer
    ssh -N -f -L localhost:8080:localhost:8080  -i [private_key] user@server.com
    ssh -N -f -L localhost:8000:localhost:8000 -i [private_key] user@server.com

```
* Navigate to the localhost port chosen
* To kill a port type in the command 
``` 
    lsof -ti:PORT_NUM | xargs kill -9
```

## Developing in Iris

Python and Javascript edits can be made in any text editor. Note any Javascript/HTML/CSS edits must be recompiled by:
```
cd [PATH_TO_IRIS]/app
webpack
```
You must restart Electron GUI interface for changes to take effect.

## Additional modularity

Our current inferential pipeline could be easily extended to a large number of other biomedical tasks, some of which may aid the other BDT teams. Additional commands can be added to iris-ncats within any text editor or the interface GUI as such:
```python
from iris import state_types as t
from iris import IrisCommand

class TreatDisease(IrisCommand):
    # what iris will call the command + how it will appear in a hint
    title = "how do I treat {disease}?"
    
    # give an example for iris to recognize the command
    examples = ["what is the treatment for {disease}", "how is {disease} treated"]
    
    # type annotations for each command argument, to help Iris collect missing values from a user
    argument_types = {"disease":t.String("Please enter the disease:")}
    
    # core logic of the command
    def command(self, disease):
        ###### CALL BACKEND ANALYSIS PIPELINE #####
        return "Answer"
        
    # wrap the output of a command to display to user
    # by default this will be an identity function
    # each element of the list defines a separate chat bubble
    def explanation(self, result):
        return ["Here is the answer", result]
```
## NCATS Specific Data and Functions

To setup the data folder: 

1. Get data as below:
```
wget 'https://www.dropbox.com/s/iwl58rre8qcav6d/DB_data.tar.gz?dl=1' -O 'DB_data.tar.gz'

```

2. Place `DB_data` folder in the folder `backend\app\user_functions\ncats`


NCATS specific functions can be found in `backend\app\user_functions\`

Things we need to do:
* Documentation of Databases (smartAPI based)
* Documentation of biomedical algorithms (logic)
* Connecting Iris to SMART API
* Streamline directory structure

## Known Issues

The following issues are recognized and are being resolved
* OSX display issues - if you have a mouse USB plugged in while running `npm start`, display proportions will be slighlty misconfigured. Try unpluggin the mouse USB then restarting Iris.
* Windows display issues - The display characteristics if using Windows platform is off, we recognize this problem and are working to fix it

Please contact us if you have any questions or issues you wish to bring to our attention!


