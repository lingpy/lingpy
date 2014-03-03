var dolgo = {};
dolgo['char dolgo_ERROR'] = 'white';
dolgo['char dolgo_V'] = '#c86464';
dolgo['char dolgo_K'] = '#c89664';
dolgo['char dolgo_P'] = '#c8c864';
dolgo['char dolgo_H'] = '#96c864';
dolgo['char dolgo_J'] = '#64c864';
dolgo['char dolgo_M'] = '#64c896';
dolgo['char dolgo_N'] = '#64c8c8';
dolgo['char dolgo_S'] = '#6496c8';
dolgo['char dolgo_R'] = '#6464c8';
dolgo['char dolgo_T'] = '#9664c8';
dolgo['char dolgo_W'] = '#c864c8';
dolgo['char dolgo_TONE'] = '#c86496';
dolgo['char dolgo_X']='#bbbbbb';
dolgo['char dolgo_GAP'] = '#bbbbbb';

var myprevious = ['','',''];

function show(what){
  var words = what.split('.'); 
  var ar = document.getElementsByClassName('taxon');
  var k=0;
  for(var i=0;i<ar.length;i++)
  {
    var taxon = ar[i].getAttribute('taxon');
    if(taxon==words[0]){
      var chars = ar[i].children[3];
      chars = chars.children[0];
      chars = chars.children[0];
      chars = chars.children[0].children;
      for(j=0;j<chars.length;j++)
      {
        var this_char = chars[j].getAttribute('num');
        if(this_char==what)
        {
          chars[j].style.backgroundColor='white';
          chars[j].style.color='red';
          var this_c = chars[j].getAttribute('class');
          k++;
        }
        else
        {
          c = dolgo[chars[j].getAttribute('class')];
          chars[j].style.backgroundColor=c;
          chars[j].style.color='white';
        }
      }
    }
    else if(taxon==myprevious[0])
    {
      var chars = ar[i].children[3];
      chars = chars.children[0];
      chars = chars.children[0];
      chars = chars.children[0].children;
      for(j=0;j<chars.length;j++)
      {
        var this_char = chars[j].getAttribute('num');
        if(this_char==myprevious[3])
        {
          chars[j].style.backgroundColor=dolgo[chars[j].getAttribute('class')];
          chars[j].style.color='white';
        }
      }
    }
  }
  var out = '<table class="popup"><tr class="popup"><th class="popup">';  
  out = out + words[0]+' <span class="'+this_c+'">'+words[1]+'</span>, context <b>'+words[2]+'</b>, '+k+' occurrences: </th><th class="popup">All occurrences of <span class="'+this_c+'">'+words[1]+'</span>: </th></tr><tr><td class="popup">';
  out = out+'<table class="display"><tr class="display"><th class="display">Language</th><th class="display">Sound</th><th class="display">Context</th><th class="display">Matches</th><th class="display">Frequency</th></tr>';
  out = out+myjson[what][0]+'</table></td><td class="popup">';
  
  // try to find the best concepts
  var concepts = myjson[what][1];
  var cstring = '<ol class="display">';
  for(var i=0;i<concepts.length;i++)
  {
  cstring = cstring+'<li class="display"><a class="display" href="#';
  cstring = cstring+concepts[i]+'">'+concepts[i]+'</a></li>';
  }
  cstring = cstring + '</ol>';
  //out = out + '<br>All occurrences of the sound:<br>';
  out = out + cstring;
  out = out + '</td></tr></table>';

  var mod = document.getElementById('display');
  mod.style.display = 'block';
  mod.innerHTML = out; 
  myprevious = [words[0], words[1], words[2], what];
}

function color(mode, what){
  // switch off general display on side-panel
  var mod = document.getElementById("display");
  mod.style.display = 'none';

  var ar = document.getElementsByClassName('taxon');
  for(var i=0;i<ar.length;i++)
  {
    var taxon = ar[i].getAttribute('taxon');
    if(taxon==myprevious[0])
    {
      var chars = ar[i].children[3];
      chars = chars.children[0];
      chars = chars.children[0];
      chars = chars.children[0].children;
      for(j=0;j<chars.length;j++)
      {
        var this_char = chars[j].getAttribute('num');
        if(this_char==myprevious[3])
        {
          chars[j].style.backgroundColor=dolgo[chars[j].getAttribute('class')];
          chars[j].style.color='white';
        }
      }
    }
  }

  // choose relevant block
  if(mode=='confidence')
  {
    var trs = document.getElementById(what).children[0].children;
    for(var i=2;i<trs.length;i++)
    {
       var chars = trs[i].children[3];
       chars = chars.children[0];
       chars = chars.children[0];
       chars = chars.children[0].children;
       if(chars[0].innerHTML != '--')
       {
         for(var k=0;k<chars.length; k++)
         {
           // get the confidence color
           var cf = chars[k].getAttribute('confidence');
           var cfc = Math.round(255-(255*cf/100));
           cfc = 'rgb('+cfc+','+cfc+','+cfc+')';
           chars[k].style.backgroundColor=cfc;
           var c = chars[k].getAttribute('char');
           chars[k].innerHTML = c;
           if(cf >= 50)
           {
             chars[k].style.color = "white";
           }
           else
           {
             chars[k].style.color = "black";
           }  
         }
       }
    }
  }
  // choose relevant block
  else if(mode=='sound_classes')
  {
    var trs = document.getElementById(what).children[0].children;
    for(var i=2;i<trs.length;i++)
    {
      var chars = trs[i].children[3];
      chars = chars.children[0];
      chars = chars.children[0];
      chars = chars.children[0].children;
      if(chars[0].innerHTML != '--' && typeof(chars) != 'undefined')
      {
        for(var k=0;k<chars.length; k++)
        {
          // get the sound class color
          var col = dolgo[chars[k].getAttribute('class')];
          var c = chars[k].getAttribute('char');
          chars[k].style.backgroundColor=col;
          chars[k].style.color = 'white';
          chars[k].innerHTML = c;
        }
      }
    }
  }
  // choose relevant block
  else if(mode=='scores')
  {
    var trs = document.getElementById(what).children[0].children;
    for(var i=2;i<trs.length;i++)
    {
      var chars = trs[i].children[3];
      chars = chars.children[0];
      chars = chars.children[0];
      chars = chars.children[0].children;
      if(chars[0].innerHTML != '--')
      {
        for(var k=0;k<chars.length; k++)
        {
          // get the confidence color
          var cf = chars[k].getAttribute('confidence');
          var cfc = Math.round(255-(255*cf/100));
          cfc = 'rgb('+cfc+','+cfc+','+cfc+')';
          chars[k].style.backgroundColor=cfc;
          chars[k].innerHTML = cf;
          if(cf >= 50)
          {
            chars[k].style.color = "white";
          }
          else
          {
            chars[k].style.color = "black";
          }  
        }
      }
    }
  }
}
