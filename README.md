# iris-ncats: A Conversational Agent for Biological Data Science 

This repository contains the prototype for our Stanford Effort Biomedical Data Translater that is able to answer the following questions:

* Given a condition, are any genetic conditions potentially protective and why?
* Given a drug and a condition, construct a pathway description that explains how the drug effects its action. Additionally, generate hypotheses of other therapeutic uses for the given drug.

This prototype is based on the IRIS, which is an open source conversational agent, originally to facilitate tasks in data science for non-programmers.  We tailored it to a biomedical context, allowing it to recognize key entities (drugs, genes, diseases), pull from various biomedical data sources, construct queries based on free text inputted by the user, and perform various analyses to produce the final answer or dossier.

The original IRIS system was created by the Bernstein lab at Stanford Computer Science, and the [iris-agents](https://github.com/Ejhfast/iris-agent) repository can be found here.

![interface](/interface_ncats.png)


## Getting Starting 

### Note: This is currently a prototype release. Additions and bux fixes are being performed.

The backend of Iris is written in Python while the front end is written in Javascript/CSS/HTML.  These are instructions to install and run iris-ncats in debugging mode that can be used for developers. We recommend using OSX for running the debugging mode. Windows and Linux versions are still in pre-release. A self-contained Electron app for OSX will be released in the next 10 months.


### Prerequisites

* [anaconda 3](https://conda.io/docs/install/quick.html). (Make sure to run `source ~/.bash_profile` after you have installed Anaconda, if it is not appearing in your path.)
* [python 2.7](https://www.python.org/download/releases/2.7/). This is in order to install the javascript components
* [xcode](https://developer.apple.com/xcode/downloads/)
* node.js: 
    * run `ruby -e "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/install)"`
    * `brew install node` (if you hit errors: go [here](http://blog.teamtreehouse.com/install-node-js-npm-mac))

### DATA (INSERT STUFF HERE!!!)

### Install and run the Python components:
```
    cd iris-agent/backend
    # create new conda env called iris with necessary packages
    conda create --name iris --file packages.txt
    # enter the conda env
    source activate iris
    # install the remaining pip packages
    pip install -r requirements.txt
```

### Install and run the Javascript components:
```
    cd iris-agent
    # If you do not already have webpack
    npm install webpack -g
    # Now installing iris javascript components
    npm install --python=python2.7
    # build JS app with webpack (can also run webpack --watch in seperate command window)
    webpack
```

## Using Iris (dev mode)
In order to run Iris in development mode. 

1. Open two terminals, one that runs frontend launcher commands, the other your Python-based backend.
2. In the python terminal type:
```
    source activate iris 


    # run the backend application
    cd iris-agent/backend/app # enter backend path

    python app.py
```
3. Wait until your backend has finished booting up (should take between 30-40 seconds)
4. In the other window, type:
```

    cd iris-agent/backend

    # start electron (this will open the application automatically)
    npm start
```

## Developing in Iris

Python and Javascript edits can be made in any text editor. Note any Javascript/HTML/CSS edits must be recompiled by:
```
cd iris-agent/app
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

## Known Issues

The following issues are recognized and are being resolved
* OSX display issues - if you have a mouse USB plugged in while running `npm start`, display proportions will be slighlty misconfigured. Try unpluggin the mouse USB then restarting Iris.
* Windows display issues - The display characteristics if using Windows platform is off, we recognize this problem and are working to fix it


