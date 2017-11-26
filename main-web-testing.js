// const electron = require('electron');
// const app = electron.app;
// const BrowserWindow = electron.BrowserWindow;



const express = require('express');
const app = express();
var https = require('https');
var http = require('http');
var path = require('path');

app.engine('html', require('ejs').renderFile);
app.set('view engine', 'html');
app.set('views', __dirname + '/app');

app.set('app', __dirname + '/app');

 
app.use(express.static(path.join(__dirname, 'app/styles')));
app.get('/', function(req, res){
  res.render('index');
});

// const port = process.env.PORT || 8000;
// const baseUrl = `http://localhost:${port}`;

// app.get('/',function(req,res){
//   const path = require('path').join(__dirname, '');

//   res.sendFile(path+'/dist/x.html');
//   //It will find and locate index.html from View or Scripts
// });

// app.on('ready', () => {
//   mainWindow = new BrowserWindow({
//     height: 800,
//     width: 1200
//   })

//   // load the local HTML file
//   let url = require('url').format({
//     protocol: 'file',
//     slashes: true,
//     pathname: require('path').join(__dirname, '/dist/index.html')
//   })
//   console.log(url)
//   mainWindow.loadURL(url)
// })

pyProc = null;

const createPyProc = () => {
  const path = require('path').join(__dirname, '');
  console.log(pyProc);
  if (pyProc === null){
    pyProc = require('child_process').spawn('bash', [path+'/'+'start_backend.sh', path])
    console.log(pyProc);
  }
  if (pyProc != null) {
    console.log('child process success')
  }
}

const exitPyProc = () => {
  console.log( "Closed out remaining connections.");
  const path = require('path').join(__dirname, '');
  killit = require('child_process').spawn('bash', [path+'/'+'stop_backend.sh'])
  // pyProc.kill()
}

// app.on('ready', createPyProc)
// app.on('will-quit', exitPyProc)
app.get(createPyProc)


var server = app.listen(3000);
console.log('Listening at localhost:3000');


//cleanup

process.on( 'SIGTERM', function () {

   server.close(exitPyProc);

   setTimeout( function () {
     console.error("Could not close connections in time, forcefully shutting down");
     process.exit(1); 
   }, 30*1000);
 });

// });
// // httpServer = http.createServer(app);
// httpServer.listen(8000,function() {
//   console.log('Port closed:)');
//   httpServer.close( exitPyProc);
// });



