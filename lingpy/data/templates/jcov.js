/* Correspondence Viewer for LingPy Cognate Judgments
 *
 * author   : Johann-Mattis List
 * email    : mattis.list@lingulist.de
 * created  : 2014-08-20 11:51
 * modified : 2014-08-24 18:43
 *
 */

/* we store namespaces of this app as JCOV */
var JCOV = {};

/* store sorting options */
JCOV.SORTED = {};
JCOV.selections = [];
JCOV.concepts = [];

JCOV.settings = {};
JCOV.settings['correspondences'] = true;
JCOV.settings['cognates'] = true;
JCOV.settings['sound_selection'] = ["K", "T", "P", "S", "N", "M", "R", "J", "W", "H", "V", "1"];
JCOV.settings['context_selection'] = ['C','c','V'];
JCOV.settings['sound_order'] = ["K", "T", "P", "S", "N", "M", "R", "J", "W", "H", "V", "1"]; /* just used for constant ordering */
JCOV.settings['refresh'] = false;
JCOV.settings['hide_empty_corrs'] = true;
JCOV.settings['collapse_cognates'] = false;
JCOV.settings['timeout'] = 10;
JCOV.settings['reduce_alignments'] = true;

JCOV.toggleRefresh = function()
{
  var tmp = document.getElementById('toggle_refresh');
  if(JCOV.settings['refresh'])
  {
    JCOV.settings['refresh'] = false;
    tmp.innerHTML = tmp.innerHTML.replace('glyphicon-ok','glyphicon-remove');
  }
  else
  {
    JCOV.settings['refresh'] = true;
    tmp.innerHTML = tmp.innerHTML.replace('glyphicon-remove','glyphicon-ok');
  }
}

JCOV.togglePanel = function(elm)
{
  $('#'+elm).toggleClass('invisible');
  var tmp = document.getElementById(elm+'_toggle');
  if(tmp.innerHTML.indexOf('glyphicon-ok') != -1)
  {
    tmp.innerHTML = tmp.innerHTML.replace('glyphicon-ok','glyphicon-remove');
    JCOV.settings[elm] = false;
  }
  else
  {
    tmp.innerHTML = tmp.innerHTML.replace('glyphicon-remove','glyphicon-ok');
    JCOV.settings[elm] = true;
  }
  
  if(JCOV.settings['refresh'])
  {
    JCOV.setHeight();
  }
}

JCOV.doIt = function(code)
{
  $('#popup_background').toggle();
  setTimeout(function(){
    code();
    $('#popup_background').toggle();
  },JCOV.settings.timeout);
}

JCOV.loading = function(start)
{
  if(start)
  {
    $('#popup_background').toggle();  
  }
  else
  {
    $('#popup_background').hide();
  }
}

/* save order of elements in ctab numbers */
JCOV.saveOrder = function(element,class_name)
{
  /* add sorted status of current table to sorter */
  var tmp = [];
  var lines = document.getElementsByClassName(class_name);
  if(lines.length != 0)
  {
    for(var i=0,elm;elm=lines[i];i++)
    {
      tmp.push(elm.dataset['value']);
    }
    JCOV.SORTED[element] = tmp;
  }  
  JCOV.resetSelection(); 
}

/* remove ordered elements from array */
JCOV.destroyOrder = function(element)
{
  delete JCOV.SORTED[element];
  JCOV.resetSelection();
}

JCOV.exportData = function(element)
{
  var elm = document.getElementById(element);
  var txt = '';
  for(var i=0,row;row=elm.rows[i];i++)
  {
    var this_row = [];
    for(var j=0,cell;cell=row.cells[j];j++)
    {
      this_row.push($(cell).text());
    }
    if(this_row.length > 1)
    {
      txt += this_row.join('\t').replace(/^\t/,'')+'\n';
    }
    else
    {
      txt += '# '+this_row[0]+'\n';
    }
  }
  
  //var store = document.getElementById('store');
  var blob = new Blob([txt], {type: 'text/plain;charset=utf-8'});
  saveAs(blob, element+'.tsv');
}

JCOV.showCorrs = function (this_lang,show)
{
  JCOV.settings['current_language'] = this_lang;

  /* just write the simple header notice to this part */
  var source = document.getElementById('correspondence_header');
  var header = '<span class="handle glyphicon glyphicon-move "></span> ';
  header += 'Frequent Sound Correspondences of '+this_lang+":";
  header += '<span title="toggle table" onclick="JCOV.togglePanel(\'correspondences\');" class="pointed glyphicon glyphicon-remove pull-right"></span>';
  header += '<span title="export current table" onclick="JCOV.exportData(\'correspondence_table\')" class="pointed glyphicon glyphicon-export pull-right"></span>';

  if(this_lang in JCOV.SORTED)
  {
    header += '<span title="destroy currently saved order of items" onclick="JCOV.destroyOrder(\''+this_lang+'\');" class="pointed glyphicon glyphicon-floppy-remove pull-right"></span> ';
  }
  header += '<span title="save current order of items" onclick="JCOV.saveOrder(\''+this_lang+'\',\'ctab_numbers\');" class="pointed glyphicon glyphicon-floppy-save pull-right"></span> ';
  source.innerHTML = header;

  var minoccs = parseInt(document.getElementById('minoccs').value);
  
  /* check for sorted stuff */
  if(this_lang in JCOV.SORTED)
  {
    var keys = JCOV.SORTED[this_lang];
  }
  else
  {
    var keys = [];
  
    /* get all keys for the language */
    for(key in CORRS)
    {
      var tmp = key.split('.')[0];
      if(tmp == this_lang)
      {
        keys.push(key);
      }
    }
    keys.sort(
        function(x,y)
        {
          var charx = x.split('.');
          var chary = y.split('.');
          var clx = getSoundClass(charx[1]);
          var cly = getSoundClass(chary[1]);
          if(clx == cly)
          {
            if(charx[1] == chary[1])
            {
              if(charx[2] == 'C')
              {
                return -1;
              }
              else
              {
                return 1;
              }
            }
            else
            {
              return charx[1].charCodeAt(0) - chary[1].charCodeAt(0);
            }
          }
          else
          {
            return JCOV.settings['sound_order'].indexOf(clx) - JCOV.settings['sound_order'].indexOf(cly);
          }
        })
  }
  
  /* get current selection of languages */
  var olangs = [];
  for(var i=0;i<LANGS.length;i++)
  {
    if(this_lang != LANGS[i] && JCOV.selections.indexOf(LANGS[i]) != -1)
    {
      olangs.push(LANGS[i]);
    }
  }
  olangs.sort();
  
  /* create the header */
  var header = '<thead>';
  header += '<tr><th class="sound_handle"></th><th class="source">SOUND</th><th class="source">CTXT</th><th class="source">OCC</th>';
  for(var i=0,olang;olang=olangs[i];i++)
  {
    header += '<th title="show correspondences for '+olang+'" class="targets pointed" onclick="JCOV.showCorrs(\''+olang+'\')">'+olang+'</th>';
  }
  header += '</tr></thead>';
  
  var txt = '<table id="correspondence_table">'+header;
  txt += '<tbody>';
  
  var counter = 0;

  for(var i=0;i<keys.length;i++)
  {
    key = keys[i];
    
    var this_sound = key.split('.')[1];
    var this_context = key.split('.')[2];
    
    var this_sound_class = getSoundClass(this_sound);
    
    if(JCOV.settings['sound_selection'].indexOf(this_sound_class) != -1 
        && JCOV.settings['context_selection'].indexOf(this_context) != -1)
    {

      var tbl = [];

      for(var j=0;j<olangs.length;j++)
      {
        tbl.push('-');
      }
      
      /* set the tracker to guarantee that empty tables are not displayed */
      var track = 0;
      for(k in CORRS[key])
      {
        var tmp = k.split('.');
        var other_lang = tmp[0];
        var other_sound = tmp[1];

        var idx = olangs.indexOf(other_lang);
        
        /* make sure sound occurs more than just once in correspondence relation */
        if(CORRS[key][k] >= minoccs)
        {
          if(tbl[idx] == '-')
          {
            tbl[idx] = plotWord(other_sound) +'<span style="display:table-cell">('+CORRS[key][k]+')</span>';
          }
          else
          {
            tbl[idx] += '<div style="margin-bottom:2px;">'+plotWord(other_sound) + '<span style="display:table-cell">('+CORRS[key][k]+')</span>';
          }
        }
        else
        {
          track += 1;
        }
      }
      if(tbl.join('').replace(/-/g,'') == '' && JCOV.settings['hide_empty_corrs']){}
      else
      {
        counter += 1;
        txt += '<tr class="correspondence_row" id="tr_'+key+'"><td data-value="'+key+'" class="ctab_numbers sound_handle pointed" id="cnum_'+counter+'"><span class="glyphicon glyphicon-move"></span></td><td title="show all occurrences of '+this_sound+'/'+this_context+'" class="source pointed" onclick="JCOV.showOccurrences(\''+key+'\',\'show\')">'+plotWord(this_sound)+'</td><td class="source">'+this_context+'</td><td class="source">'+OCCS[key].length+'</td><td>'+tbl.join('</td><td>')+'</td></tr>';
      }
    }
  }
  
  /* make table sortable */
  document.getElementById('correspondence_data').innerHTML = txt+'</tbody></table>';
  
  $('#correspondence_table tbody').sortable({
    start: function( event, ui ) {      
      clone = $(ui.item[0].outerHTML).clone();},
    placeholder: {
      element: function(clone, ui) {
        return $('<tr class="selected" style="opacity:0.2;">'+clone[0].innerHTML+'</tr>');},
      update: function() {
        return;}},
    handle: ".ctab_numbers"
  });
  
  /* initialize correspondence table as fixed-header-table */
  $('#correspondence_table').fixedHeaderTable({ 
  });

  if(typeof show != 'undefined')
  {
    if(document.getElementById('correspondences').className.indexOf('invisible') != -1)
    {
      JCOV.togglePanel('correspondences');
      JCOV.setHeight();
    }
  }
}

JCOV.collapseCognates = function (elm)
{
  var trs = document.getElementsByClassName(elm.id);
  if(elm.dataset['collapsed'] == 'false')
  {
    /* note that we store the elements width here in order to
     * guarantee that the element's width won't change too much 
     * if all items in the table are collapsed 
     */
    var width = $(elm).width(); //offsetWidth;
    elm.dataset['collapsed'] = true;
    for(var i=0,tr;tr=trs[i];i++)
    {
      $(tr).addClass('invisible');
    }
    elm.innerHTML = elm.innerHTML.replace('collapse-down','collapse-up');
    var nwidth = $(elm).width();
    if(width > nwidth+10)
    {
      $(elm).width(width);
    }
  }
  else
  {
    elm.dataset['collapsed'] = 'false';
    for(var i=0,tr;tr=trs[i];i++)
    {
      $(tr).removeClass('invisible');
    }
    elm.innerHTML = elm.innerHTML.replace('collapse-up','collapse-down');
  }
}

/* function displays all occurrences of a given word, or alternatively,
 * all cognate sets in the COGS-table */
JCOV.showOccurrences = function (sound,show)
{
  /* set up the header */
  var header = '<span class="glyphicon glyphicon-move handle"></span> ';
  var occs;

  /* store as current selection */
  JCOV.settings['current_sound'] = sound;

  /* check for undefined occurrences */
  if(sound == 'ALL')
  {
    var clen = Object.keys(WLS).length;

    /* check whether concepts are of zero length */
    if('ALL' in JCOV.SORTED)
    {
      occs = JCOV.SORTED['ALL'];
    }
    else if(JCOV.concepts.length == 0)
    {
      occs = Object.keys(WLS);
      occs.sort();
      occs = occs.slice(0,20);
      JCOV.concepts = occs;
    }
    else
    {
      occs = JCOV.concepts;
    }
    header += 'Showing cognate sets occurring in '+occs.length + ' out of '+clen+' concepts:';
  }
  else
  {
    if(sound in JCOV.SORTED)
    {
      occs = JCOV.SORTED[sound];
    }
    else
    {
      occs = OCCS[sound];
    }

    /* get all three values*/
    var tmp = sound.split('.');
    var tlang = tmp[0];
    var tsound = tmp[1];
    var tcontext = tmp[2];
    
    if(occs.length == 1)
    {
      header += tlang+' '+plotWord(sound.split('.')[1],'div')+ '/'+tcontext+' occurs in '+occs.length+' concept:';
    }
    else
    {
      header += tlang+' '+plotWord(sound.split('.')[1],'div')+ '/'+tcontext+' occurs in '+occs.length+' concepts:';
    }
  }
  /* fill in the rest of the header */
  header += '<span title="hide panel" onclick="JCOV.loading(true);setTimeout(function(){JCOV.togglePanel(\'cognates\');},JCOV.settings.timeout);JCOV.loading(false);" class="glyphicon glyphicon-remove pull-right pointed"></span>';
  header += '<span title="show all cognate sets" onclick="JCOV.doIt(function(){JCOV.showOccurrences(\'ALL\');});" class="glyphicon glyphicon-list pull-right pointed"></span>';
  if(JCOV.settings['collapse_cognates'])
  {
    header += '<span title="expand all" onclick="JCOV.doIt(function(){JCOV.settings[\'collapse_cognates\']=false;JCOV.showOccurrences(\''+sound+'\',\''+show+'\');});" class="glyphicon glyphicon-resize-full pointed pull-right"></span>';
  }
  else{
    header += '<span title="collapse all" onclick="JCOV.doIt(function(){JCOV.settings[\'collapse_cognates\']=true;JCOV.showOccurrences(\''+sound+'\',\''+show+'\');});" class="glyphicon glyphicon-resize-small pull-right pointed"></span>';
  }

  header += '<span title="export current table" onclick="JCOV.exportData(\'cognates_table\')" class="pointed glyphicon glyphicon-export pull-right"></span>';

  if(sound in JCOV.SORTED)
  {
    header += '<span title="destroy currently saved order of items" onclick="JCOV.doIt(function(){JCOV.destroyOrder(\''+sound+'\');});" class="pointed glyphicon glyphicon-floppy-remove pull-right"></span> ';
  }
  header += '<span title="save current order of items" onclick="JCOV.doIt(function(){JCOV.saveOrder(\''+sound+'\',\'cognates_table_header\');});" class="pointed glyphicon glyphicon-floppy-save pull-right"></span> ';

  var txt = '<table id="cognates_table">';
  txt += '<thead><tr><th>ID</th><th>Doculect</th><th>Concept</th><th>IPA</th><th>COGID</th><th>Alignment</th></tr></thead>';
  
  for(var i=0,concept;concept=occs[i];i++)
  {
    txt += '<tbody><tr><th id="gloss_'+GlossId[concept]+'" data-value="'+concept+'" ';
    if(JCOV.settings['collapse_cognates'])
    {
      txt += 'data-collapsed="true"';
      txt += ' class="cognates_table_header" colspan="6">&quot;'+concept+'&quot;';
      txt += '<span class="pointed glyphicon glyphicon-collapse-down pull-right" onclick="JCOV.collapseCognates(this.parentNode)"></span>';
      txt += '<span class="pointed handlerx glyphicon glyphicon-move pull-left"></span> '; 
    }
    else
    {
      txt += 'data-collapsed="false"';
      txt += ' class="cognates_table_header" colspan="6">&quot;'+concept+'&quot;';
      txt += '<span class="pointed glyphicon glyphicon-collapse-up pull-right" onclick="JCOV.collapseCognates(this.parentNode)"></span>';
      txt += '<span class="pointed handlerx glyphicon glyphicon-move pull-left"></span> '; 
    }

    txt += '</th></tr>';
    
    /* quick search for cognate words */
    var entries = WLS[concept];
    
    var gcid = '';
    for(var j=0,entry;entry=entries[j];j++)
    {
      if(entry[0] == tlang)
      {
        if(entry[2].indexOf(tsound) != -1)
        {
          gcid = entry[3];
        }
      }
    }
    
    var pcogid = '';
    var pcol = 'white';
    var alms = [];
    
    var alm_entry = [];
    
    /* get entries matching language selection condition */
    var nentries = [];
    for(var j=0,entry;entry=entries[j];j++)
    {
      if(JCOV.selections.indexOf(entry[0]) != -1)
      {
        nentries.push(entry);
        alm_entry.push(entry[4]);
      }
    }
    
    if(JCOV.settings['reduce_alignments'])
    {
      /* iterate over entries and reduce alignments */
      for(var j=0,entry;entry=nentries[j];j++)
      {
        var tcogid = entry[3];
        /* first alignment is alway special */
        if(j == 0)
        {
          alms.push(entry[4].split(' '));
          pcogid = tcogid;
        }
        /* last alignmetn needs also special treatment */
        else if(j == (nentries.length-1))
        {
          /* if last alignment is different from alignments before, first insert
           * preceding alignments, then insert last alignment */
          if(tcogid != pcogid)
          {
            alms = transpose(alms);
            var new_alms = [];
            for(var k=0;k<alms.length;k++)
            {
              if(alms[k].join('').replace(/-/g,'') == '')
              {}
              else
              {
                new_alms.push(alms[k]);
              }
            }

            alms = transpose(new_alms);

            for(var k=0;k<alms.length;k++)
            {
              alm_entry[(j-(alms.length-k))] = alms[k].join(' ');
            }

            nentries[j][4] = entry[4].replace(/\s*-/g,' ');
          }
          else
          {
            alms.push(entry[4].split(' '));
            alms = transpose(alms);
            var new_alms = [];
            for(var k=0;k<alms.length;k++)
            {
              if(alms[k].join('').replace(/-/g,'') == '')
              {}
              else
              {
                new_alms.push(alms[k]);
              }
            }

            alms = transpose(new_alms);

            for(var k=0;k<alms.length;k++)
            {
              alm_entry[(j-(alms.length-k)+1)] = alms[k].join(' ');
            }
          }
        }
        else if(tcogid != pcogid)
        {
              alms = transpose(alms);
              var new_alms = [];
              for(var k=0;k<alms.length;k++)
              {
                if(alms[k].join('').replace(/-/g,'') == '')
                {
                }
                else
                {
                  new_alms.push(alms[k]);
                }
              }
              alms = transpose(new_alms);
              for(var k=0;k<alms.length;k++)
              {
                alm_entry[(j-(alms.length-k))] = alms[k].join(' ');
              }
              alms = [];
              alms.push(entry[4].split(' '));
          pcogid = tcogid;
        }
        else
        {
          alms.push(entry[4].split(' '));
        }
      }
    }

    for(var j=0,entry;entry=nentries[j];j++)
    {
      
      var tcogid = entry[3];

        /* check for previously differing cogid */
        if(tcogid != pcogid)
        {
          pcogid = tcogid;
          if(pcol == 'white')
          {
            pcol = 'lightgray';
          }
          else
          {
            pcol = 'white';
          }
        }
        
        /* collapse cognates check */
        if(JCOV.settings['collapse_cognates'])
        {
          txt += '<tr class="invisible gloss_';
        }
        else
        {
          txt += '<tr class="gloss_';
        }
        txt += GlossId[concept]+'" ';

        if(gcid == tcogid)
        {
          txt += 'style="background-color:LightBlue">';
        }
        else
        {
          txt += 'style="background-color:'+pcol+'">';
        }
        txt += '<td>'+entry[1]+'</td>'; // id
        txt += '<td style="cursor:pointer" title="show sound correspondences for this doculect" onclick="JCOV.showCorrs(\''+entry[0]+'\',\'show\');">'+entry[0]+'</td>'; // language
        txt += '<td>'+concept+'</td>'; // concept
        txt += '<td>'+entry[2]+'</td>'; // ipa
        txt += '<td>'+entry[3]+'</td>'; // cogid
        txt += '<td>'+plotWord(alm_entry[j])+'</td>'; // alignment
        txt += '</tr>';
      //}
    }
    txt += '</tbody>';
  }
  
  txt += '</table>';
  
  document.getElementById('cognates_header').innerHTML = header;
  document.getElementById('cognates_data').innerHTML = txt;

  $('#cognates_table').sortable({
    group: '.blu',
    handle: '.handlerx',
    start: function( event, ui ) {      
      clone = $(ui.item[0].outerHTML).clone();},
    placeholder: {
      element: function(clone, ui) {
        return $('<tbody class="selected" style="opacity:0.2;">'+clone[0].innerHTML+'</tbody>');},
    update: function() {
      return;}}
   });
  
  $('.sortables').addClass('ui-helper-clearfix');
  
  /* initialize correspondence table as fixed-header-table */
  $('#cognates_table').fixedHeaderTable({});
  JCOV.setHeight();

  if(typeof show != 'undefined')
  {
    if(document.getElementById('cognates').className.indexOf('invisible') != -1)
    {
      JCOV.togglePanel('cognates');
      JCOV.setHeight();
    }
    else
    {
      JCOV.setHeight();
    }
  }

}

JCOV.createSelectors = function()
{
  /* create concept selector */
  var tmp_selections = [];
  var txt = '<select class="multiselect" id="selectora" multiple>';
  var counter = 1;
  var concepts = Object.keys(WLS);
  concepts.sort();

  for(var i=0,concept;concept=concepts[i];i++)
  {
    if(counter < 21)
    {
      tmp_selections.push(concept);
      txt += '<option value="'+concept+'" selected>'+concept+'</option>';
    }
    else
    {
      txt += '<option value="'+concept+'">'+concept+'</option>';
    }
    counter += 1;
  }
  txt += '</select>';
  document.getElementById('concept_selector').innerHTML = txt;

  /* create language selector */
  var tmp_selections = [];
  var txt = '<select class="multiselect" id="selectorix" multiple>';
  var counter = 1;
  for(var i=0,lang;lang=LANGS[i];i++)
  {
    if(counter < 7)
    {
      tmp_selections.push(lang);
      txt += '<option value="'+lang+'" selected>'+lang+'</option>';
    }
    else
    {
      txt += '<option value="'+lang+'">'+lang+'</option>';
    }
    counter += 1;

  }
  txt += '</select>';
  document.getElementById('language_selector').innerHTML = txt;
  
  $('#selectorix').multiselect({
    selectAll:true,
    enableFiltering:true,
    maxHeight: window.innerHeight-150, /* change later to cog_height */
    buttonClass: 'btn-link',
    includeSelectAllOption: true,
    enableCaseInsensitiveFiltering: true,
    buttonContainer: '<div style="display:inline" />',
    buttonText: function(options,select){return 'Select Doculects <b class="caret"></b>';}
  }); 
  $('#selectora').multiselect({
    selectAll:true,
    enableFiltering:true,
    maxHeight: window.innerHeight-150, /* change later to cog_height */
    buttonClass: 'btn-link',
    includeSelectAllOption: true,
    enableCaseInsensitiveFiltering: true,
    buttonContainer: '<div style="display:inline" />',
    buttonText: function(options,select){return 'Select Concepts <b class="caret"></b>';}
  }); 

  
  JCOV.selections = tmp_selections;

  /* create sound selector */
  var usym_map= {
    "K":"Velars, palatal affricates",
    "1":"Tones",
    "H":"Glottals and laryngeals",
    "T":"Dentals",
    "V":"Vowels",
    "N":"Nasals",
    "P":"Labials",
    "R":"Liquids",
    "S":"Sibilants",
    "W":"Labial glide",
    "M":"Labial nasal",
    "J":"Palatal glide",
  };
  var usyms = ["P", "T", "K", "M", "N", "R", "S", "J", "W", "H", "V", "1"];
  /* unique symbols in DOLGO */
  
  /* make select button */
  var txt = '<select class="multiselect" id="selectorway" multiple="multiple">';
  txt += '<optgroup label="Contexts">';
  txt += '<option class="context" value="/C" selected>Sonority increases</option>';
  txt += '<option class="context" value="/c" selected>Sonority decreases</option>';
  txt += '<option class="context" value="/V" selected>Sonority peak</option>';
  txt += '</optgroup>';
  txt += '<optgroup label="Sound Classes">';
  
  for(var i=0,sym;sym=usyms[i];i++)
  {
    txt += '<option value="'+sym+'" selected>'+usym_map[sym]+'</option>';
  }
  txt += '</option>';
  txt += '</optgroup></select>';

  document.getElementById('sound_selector').innerHTML = txt;

  $('#selectorway').multiselect({
    selectAll:true,
    enableFiltering:true,
    maxHeight: window.innerHeight-150, /* change later to cog_height */
    buttonClass: 'btn-link',
    label: function(elm){if(elm.className == "context"){return elm.innerHTML;} else{return '<div class="residue dolgo_'+elm.value+'" style="font-weight:normal;padding:2px;font-family:Sans-Serif;">'+elm.innerHTML+'</div>';}},
    includeSelectAllOption: true,
    enableCaseInsensitiveFiltering: true,
      buttonContainer: '<div style="display:inline" />',
      buttonText: function(options,select){return 'Select Sounds <b class="caret"></b>';}
  }); 
  
  $('#selectorz').multiselect({
    selectAll:true,
    maxHeight: window.innerHeight-150, /* change later to cog_height */
    buttonClass: 'btn-link',
    label: function(elm){return elm.innerHTML;},
    includeSelectAllOption: true,
    buttonContainer: '<div style="display:inline" />',
    buttonText: function(options,select){return 'Settings <b class="caret"></b>';},
    onChange: function(elm, checked) {if(checked) { JCOV.settings[elm[0].value] = true; } else { JCOV.settings[elm[0].value] = false; } } });

  $('.multiselect-container')
    .css('overflow-y','auto')
    .css('padding',"15px");
}

JCOV.resetSelection = function()
{
  /* reset concepts selection */
  var selector = document.getElementById('selectora');
  var concept_selections = [];
  for(var i=0,option;option=selector.options[i];i++)
  {
    if(option.selected)
    {
      concept_selections.push(option.value);
    }
  }
  JCOV.concepts = concept_selections;
  
  /* reset sounds selection */
  var selector = document.getElementById('selectorway');
  var sound_selections = [];
  var context_selections = [];
  for(var i=0,option;option=selector.options[i];i++)
  {
    if(option.selected && option.value.slice(0,1) == '/')
    {
      context_selections.push(option.value.slice(1,2));
    }
    else if(option.selected)
    {
      sound_selections.push(option.value);
    }
  }
  if(sound_selections.length > 0)
  {
    JCOV.settings['sound_selection'] = sound_selections;
  }
  if(context_selections.length > 0)
  {
    JCOV.settings['context_selection'] = context_selections;
  }
  
  /* reset languages selection */
  var selector = document.getElementById('selectorix');
  var new_selections = [];
  for(var i=0,option;option=selector.options[i];i++)
  {
    if(option.selected)
    {
      new_selections.push(option.value);
    }
  }

  /* check for currently visible panels */
  var corrs = document.getElementById('correspondences');
  var cogs = document.getElementById('cognates');
  if(cogs.className.indexOf('invisible') != -1)
  {
    cogs = false;
  }
  else
  {
    cogs = true;
  }
  if(corrs.className.indexOf('invisible') != -1)
  {
    corrs = false;
  }
  else
  {
    corrs = true;
  }
  
  /* check for validity of selection */

  if(new_selections.indexOf(JCOV.settings['current_language']) != -1)
  {
    var main_lang = JCOV.settings['current_language'];
  }
  else if(new_selections.length > 0)
  {
    var main_lang = new_selections[0];
  }
  else
  {
    new_selections = LANGS;
    var main_lang = new_selections[0];
  }

  /* store selections */
  JCOV.selections = new_selections;

  if(cogs && corrs)
  {
    JCOV.showCorrs(main_lang);
    JCOV.showOccurrences(JCOV.settings['current_sound']);
  }
  else if(cogs)
  {
    if(JCOV.settings['current_sound'])
    {
      JCOV.showOccurrences(JCOV.settings['current_sound']);
    }
    else
    {
      JCOV.showOccurrences('ALL');
    }
  }
  else if(corrs)
  {
    JCOV.showCorrs(main_lang);
  }
  
  JCOV.setHeight();
}

JCOV.setHeight = function()
{
  var cur_height = window.innerHeight;
  var cur_width = window.innerWidth;

  var cor_height = cur_height - 140;
  
  if(JCOV.settings['correspondences'] && JCOV.settings['cognates'])
  {
    var cor_width = cur_width * 0.5;
    var cog_width = cur_width * 0.45;
  }
  else if(JCOV.settings['correspondences'])
  {
    var cor_width = cur_width * 0.9;
    var cog_width = cur_width * 0.45;
  }
  else
  {
    var cor_width = cur_width * 0.9;
    var cog_width = cur_width * 0.9;
  }

  $('#correspondence_data').css('max-width',cor_width+'px');
  $('#correspondence_data').css("height",cor_height+'px');
  
  
  if(JCOV.settings['correspondences'])
  {
    var cog_height = document.getElementById('correspondence_data').offsetHeight; /* check this, it sounds strange */
  }
  else
  {
    var cog_height = cor_height;
  }
  
  $('#cognates_data').css('height',cog_height+'px');
  $('#cognates_data').css('max-width',cog_width+'px');
  
  JCOV.settings['cth'] = cor_height;
  JCOV.settings['ctw'] = cur_width;

  /* resize the multi-selector */
  $('.multiselect-container').css('max-height',cur_height-100); 

  JCOV.showCorrs(JCOV.settings['current_language']);
}

JCOV.filterData = function(event,elm)
{
  if(event.keyCode != 13)
  {
    return;
  }
  else
  {
    var tmpA = elm.value.split(',');
    var tmpB = [];
    for(var i=0,val;val=tmpA[i];i++)
    {
      val = val.replace(/^\s*/,'').replace(/\s*$/,'');
      if(val in WLS)
      {
        tmpB.push(val);
      }
    }
    if(tmpB.length > 0)
    {
      JCOV.concepts = tmpB;
      JCOV.showOccurrences(JCOV.settings['current_sound']);
    }
    else
    {
      elm.value = '';
    }
  }
}

/* function resizes data upon change of window-height */
$(document).ready(function(){
  function setHeight()
  {
    JCOV.setHeight();
  }
  $(window).on('resize', function(){setHeight()});
  setHeight();
});

$('#handles').sortable({
  handle: '.handle',
  start: function( event, ui ) {      
    clone = $(ui.item[0].outerHTML).clone();},
  placeholder: {
    element: function(clone, ui) {
      return $('<div class="selected" style="opacity:0.2;">'+clone[0].innerHTML+'</div>');},
  update: function() {
    return;}}
  });

JCOV.init = function()
{

  JCOV.createSelectors();
  JCOV.showOccurrences('ALL');
  JCOV.togglePanel('correspondences');
  JCOV.showOccurrences('ALL');
  

  var stuff = Object.keys(WLS);

  function split(val ) {
    return val.split(/,\s*/);
  }
  function extractLast(term ) {
    return split(term).pop();
  }
  
  $('#filter_data').bind('keydown', function(event ) 
      {
        if (event.keyCode === $.ui.keyCode.TAB && $(this).data('ui-autocomplete').menu.unhidden)
        {
          event.preventDefault();
        }
      }
    ).autocomplete(
    {
      delay: 0,
      minLength: 0,
      source: function(request, response ) 
      {
        // delegate back to autocomplete, but extract the last term
        response($.ui.autocomplete.filter(
        stuff, extractLast(request.term)));
      },
      focus: function() 
      {
        // prevent value inserted on focus
        return false;
      },
      select: function(event, ui ) 
      {
        var terms = split(this.value);
        // remove the current input
        terms.pop();
        // add the selected item
        terms.push(ui.item.value);
        // add placeholder to get the comma-and-space at the end
        terms.push('');
        this.value = terms.join(', ');
        return false;
      }
    }
  );
  var pop = document.getElementById('popup_background');
  var spinner = new Spinner().spin();
  pop.appendChild(spinner.el);
  var brand = document.getElementById('brand');
  brand.innerHTML = 'jCoV &lt;'+FILE+'&gt;';
  brand.href = FILE;
}

JCOV.init();

$(window).load(function(){$("#popup_background").fadeOut("slow");});
