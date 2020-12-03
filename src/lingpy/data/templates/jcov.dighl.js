/* define keywords for general MSA format */
var keywords = {
  "LOCAL": "", 
  "MERGE": "", 
  "SWAPS": "", 
  "ID": "",
  "IGNORE":"",
  "COLUMNID":"",
  "COMPLEX":"",
  "STANDARD":""
};

var params = ["view", "edit", "save", "refresh", "undo_button", "redo_button"];

var privates = [
  "width",
  "uniques",
  "sequences",
  "mode",
  "taxlen",
  "type"
  ];

/* define the converter object for the coloring of the cols */
var DOLGO = {
  "t\u035c\u0255": "K", 
  "\u2082\u2085": "1", 
  "\u2082\u2084": "1", 
  "\u2082\u2083": "1", 
  "\u2082\u2082": "1", 
  "\u2082\u2081": "1", 
  "\u0127": "H", 
  "\u02a4": "K", 
  "\u1e25": "H", 
  "d\u0361\u0292": "K", 
  "d\u0361\u0291": "K", 
  "\u0221": "T", 
  "\u2083\u2082": "1", 
  "\u2085": "1", 
  "\u0277": "V", 
  "\u0235": "N", 
  "\u03b2": "P", 
  "d\u035cz": "K", 
  "\u012b": "V", 
  "\u028c": "V",
  "ʮ":"V",
  "u": "V", 
  "\u02a8": "K", 
  "\u0288": "T", 
  "t\u035cs": "K", 
  "\u1d07": "V", 
  "\u2084": "1", 
  "t\u0361\u0283": "K", 
  "t\u0361\u0282": "K", 
  "\u0280": "R", 
  "s": "S", 
  "\u011b": "V", 
  "\u0294": "H", 
  "\u2c71": "W", 
  "\u0113": "V", 
  "\u0290": "S", 
  "\u01ce": "V", 
  "\u00ec": "V", 
  "\u026d": "R", 
  "\u025e": "V", 
  "\u00e8": "V", 
  "i": "V", 
  "d": "T", 
  "\u0167": "T", 
  "\u0271": "M", 
  "e": "V", 
  "\u00e0": "V", 
  "\u0261": "K", 
  "\u027d": "R", 
  "p\u035cf": "P", 
  "\u0279": "R", 
  "\u00f4": "V", 
  "\u2075": "1", 
  "\u00f0": "T", 
  "q": "K", 
  "t\u035c\u0282": "K", 
  "\u014b": "N", 
  "\u00b3\u00b2": "1", 
  "\u00b3\u00b3": "1", 
  "\u00b3\u00b9": "1", 
  "\u0259": "V", 
  "\u2081\u2084": "1", 
  "\u2081\u2085": "1", 
  "\u2081\u2082": "1", 
  "\u2081\u2083": "1", 
  "\u2081\u2081": "1", 
  "\u026e": "R", 
  "\u1d00": "V", 
  "\u2083\u2085": "1", 
  "\u2084\u2085": "1", 
  "\u2084\u2084": "1", 
  "\u02a5": "K", 
  "\u2084\u2081": "1", 
  "\u026c": "R", 
  "\u2084\u2083": "1", 
  "\u2084\u2082": "1", 
  "\u00b9": "1", 
  "\u026a": "V", 
  "\u0234": "R", 
  "\u00b3\u2074": "1", 
  "\u00b3\u2075": "1", 
  "\u028d": "M", 
  "p\u0361f": "P", 
  "\u0289": "V", 
  "\u0285": "V", 
  "\u0281": "R", 
  "f": "P", 
  "\u029d": "S", 
  "\u0299": "W", 
  "\u0295": "H", 
  "\u0264": "V", 
  "\u025b": "V", 
  "\u0253": "P", 
  "\u2075\u00b2": "1", 
  "l": "R", 
  "\u00ed": "V", 
  "\u2083\u2081": "1", 
  "h": "H", 
  "\u00e9": "V", 
  "\u0262": "K", 
  "y": "V", 
  "\u00e1": "V", 
  "\u017e": "K", 
  "\u0278": "P", 
  "t\u035c\u0283": "K", 
  "\u0274": "N", 
  "\u00f5": "V", 
  "t\u0361s": "K", 
  "p": "P", 
  "d\u035c\u0292": "K", 
  "\u0255": "S", 
  "d\u0361z": "K", 
  "\u2074": "1", 
  "\u03c7": "K", 
  "\u00f8": "V", 
  "\u0142": "R", 
  "\u0153": "V", 
  "\u025c": "V", 
  "d\u035c\u0290": "K", 
  "d\u035c\u0291": "K", 
  "\u0258": "V", 
  "\u2083\u2084": "1", 
  "\u0254": "V", 
  "\u0251": "V", 
  "\u0152": "V", 
  "\u0250": "V", 
  "\u00b9\u2074": "1", 
  "\u00b9\u2075": "1", 
  "\u0275": "V", 
  "\u00b9\u00b9": "1", 
  "\u0129": "V", 
  "\u02a6": "K", 
  "\u2082": "1", 
  "\u0291": "S", 
  "\u2070": "1", 
  "\u0276": "V", 
  "\u03b8": "T", 
  "t\u0361\u03b8": "K", 
  "d\u0361\u0290": "K", 
  "\u2083": "1", 
  "\u00b2": "1", 
  "t": "T", 
  "\u0131": "V", 
  "\u028e": "R", 
  "\u2075\u00b3": "1", 
  "\u010d": "K", 
  "\u028a": "V", 
  "\u0272": "N", 
  "\u2075\u00b9": "1", 
  "\u0282": "S", 
  "\u0180": "P", 
  "\u0101": "V", 
  "\u0270": "J", 
  "\u0292": "S", 
  "\u0111": "T", 
  "\u00ee": "V", 
  "\u2085\u2083": "1", 
  "\u2085\u2081": "1", 
  "\u00ea": "V", 
  "k": "K", 
  "\u2085\u2084": "1", 
  "\u2085\u2085": "1", 
  "\u00e6": "V", 
  "\u0267": "S", 
  "\u00e2": "V", 
  "c": "K", 
  "\u0161": "S", 
  "\u00fe": "T", 
  "\u027f": "V", 
  "\u207b": "1", 
  "w": "W", 
  "\u014d": "V", 
  "\u00f2": "V", 
  "\u0273": "N", 
  "\u00b2\u2075": "1", 
  "\u00b2\u2074": "1", 
  "\u2074\u00b9": "1", 
  "\u2075\u2074": "1", 
  "\u2074\u00b3": "1", 
  "\u2074\u00b2": "1", 
  "\u025f": "K", 
  "_": "_", 
  "\u0257": "T", 
  "\u2083\u2083": "1", 
  "\u2081": "1", 
  "\u01d0": "V", 
  "x": "K", 
  "\u2085\u2082": "1", 
  "\u026f": "V", 
  "\u02a7": "K", 
  "\u02a3": "K", 
  "t\u0361\u0255": "K", 
  "m": "M", 
  "\u0236": "T", 
  "\u00b3": "1", 
  "o": "V", 
  "\u026b": "R", 
  "\u028f": "V", 
  "\u00b2\u00b3": "1", 
  "\u028b": "W", 
  "\u2074\u2075": "1", 
  "\u2074\u2074": "1", 
  "\u0283": "S", 
  "\u00b2\u00b9": "1", 
  "\u029f": "R", 
  "\u1d1c": "V", 
  "g": "K", 
  "\u00b9\u00b2": "1", 
  "n": "N", 
  "\u0265": "J", 
  "j": "J", 
  "\u00b9\u00b3": "1", 
  "\u0266": "H", 
  "\u00e7": "S", 
  "b": "P", 
  "\u00e3": "V", 
  "\u0263": "K", 
  "\u027e": "R", 
  "z": "S", 
  "v": "W", 
  "+": "+", 
  "a": "V", 
  "r": "R", 
  "\u00f3": "V", 
  "\u00b2\u00b2": "1", 
  "\u0148": "N", 
  "\u2075\u2075": "1", 
  "\u0144": "N", 
  "t\u035c\u03b8": "K", 
  "\u0268": "V", 
  "\u01dd": "V", 
  "\u025a": "V", 
  "\u0256": "T", 
  "\u0252": "V", 
  "\u027b": "R",
  "Ɂ": "H"
};

/* simple helper function to retrieve sound classes */
function getSoundClass(sound)
{
    
		if (sound in DOLGO){dolgo = DOLGO[sound]}
    else if (sound.slice(0,2) in DOLGO){dolgo = DOLGO[sound.slice(0,2)];}
    else if (sound.slice(0,1) in DOLGO){dolgo = DOLGO[sound.slice(0,1)];}
    else if (sound.slice(1,3) in DOLGO){dolgo = DOLGO[sound.slice(1,3)];}
    else if (sound.slice(1,2) in DOLGO){dolgo = DOLGO[sound.slice(1,2)];}
    else if (sound == "-"){dolgo = "-";}

    return dolgo;
}

function plotWord(word, tag)
{
  if(typeof tag == 'undefined')
  {
    tag = 'span';
  }
  
  try
  {
	  var phones = word.split(' ');
  }
  catch(e)
  {
    alert(word);
    return '';
  }
	var text = '';
	for(var i=0;i<phones.length;i++)// in phones)
	{
		var phon = phones[i];

    /* now try to find the column */
    var dolgo = "dolgo_ERROR";
    
		if (phon in DOLGO){dolgo = "dolgo_"+DOLGO[phon]}
    else if (phon.slice(0,2) in DOLGO){dolgo = "dolgo_"+DOLGO[phon.slice(0,2)];}
    else if (phon.slice(0,1) in DOLGO){dolgo = "dolgo_"+DOLGO[phon.slice(0,1)];}
    else if (phon.slice(1,3) in DOLGO){dolgo = "dolgo_"+DOLGO[phon.slice(1,3)];}
    else if (phon.slice(1,2) in DOLGO){dolgo = "dolgo_"+DOLGO[phon.slice(1,2)];}
    else if (phon == "-"){dolgo = "dolgo_GAP";}
    
    if(phon != '-')
    {            
	    text += '<'+tag+' class="residue '+dolgo+'">'+phon+' </'+tag+'>';
    }
		else
    {
	    text += '<'+tag+' class="residue '+dolgo+'">'+phon+' </'+tag+'>';
    }
  }

	return text;
}

function plotMorphemes(word, tag, sep)
{
  if(typeof sep == 'undefined')
  {
    sep = '\\+';
  }
  
  var text_lines = [];
  var morphemes = word.split(new RegExp('\\s'+'*'+sep+'\\s'+'*'));
  for(var i=0,m;m=morphemes[i];i++)
  {
    var morpheme = '<span class="morpheme">'+plotWord(m,tag)+'</span>';
    text_lines.push(morpheme);
  }
  return text_lines.join('<span class="boundary">.</span>');
}
function sum(list)
{
  function add(prev, foll){return prev + foll;}
  return list.reduce(add);
}
function max(list)
{
  function compare(prev,foll)
  {
    if(prev == foll){return prev;}
    else if(prev > foll){return prev;}
    else{return foll}
  }
  return list.reduce(compare);
}
function range(start,stop,step)
{
  if(typeof step == 'undefined'){step=1}
  if(typeof stop == 'undefined'){stop=start;start=0}
  
  if(stop > start && step < 0){return undefined}
  if(step == 0){return undefined}
  
  
  var list = [];
  if(start < stop)
  {
    for(var i=start;i<stop;i+=step)
    {
      list.push(i);
    }
  }
  else if(stop < start)
  {
    for(var i=start;i>stop;i+=step)
    {
      list.push(i);
    }
  }
  else{return undefined}
  return list;
}

function loadFile(url,async)
{
  if(typeof async == 'undefined')
  {
    async = false;
  }
  else
  {
    async = true;
  }

  var store = document.getElementById('store');
  if(store === null)
  {
    var store = document.createElement('div');
    store.style.display = "none";
    document.body.appendChild(store);
  }

  $.ajax({
    async: async,
    type: "GET",
    url: url,
    dataType: "text",
    success: function(data) { store.innerText = data; }
  });
}

/* http://stackoverflow.com/questions/1141302/is-there-a-sleep-function-in-javascript */
function sleep(delay) {
  var start = new Date().getTime();
  while (new Date().getTime() < start + delay);
}

function transpose(array)
{
  try{
  var new_array = [];
  for(var i=0;i<array[0].length;i++)
  {
    var tmp_row = [];
    for(var j=0;j<array.length;j++)
    {
      tmp_row.push(array[j][i]);
    }
    new_array.push(tmp_row);
  }
  return new_array;
  }
  catch(e)
  {
    alert("TRANSPOSE ERROR: "+JSON.stringify(array));
    return array;
        }
}
