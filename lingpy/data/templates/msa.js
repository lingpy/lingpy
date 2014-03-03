var sorted = "alpha";

function sortTable(){
  var tbl = document.getElementById("msa").tBodies[0];
  var store = [];
  for(var i=0, len=tbl.rows.length; i<len; i++){
    var row = tbl.rows[i];
    if(sorted=="alpha"){store.push([row.getAttribute('sequence'),row])}
    else{store.push([row.getAttribute('taxon'), row])}
  }

  store.sort(function(a,b){return a[0]-b[0];});
  
  for(var i=0, len=store.length; i<len; i++)
  {
    tbl.appendChild(store[i][1]);
  }
  store = null;
  
  if(sorted=='alpha'){sorted='sequence';}
  else{sorted='alpha';}

  var test=document.getElementById('test');
  test.innerHTML = sorted;
}
function toggle(){
  var tbl = document.getElementById("msa").tBodies[0];
  
  /* Go through all elements and set them to visible if they have not been visible before */
  for(var i=0; i<tbl.rows.length; i++)
  {
    var row = tbl.rows[i];
    if(row.getAttribute('unique') == "false")
    {
      var display = row.style.display;
      if(display != 'none')
      {
        tbl.rows[i].style.display = 'none';
      }
      else
      {
        tbl.rows[i].style.display = 'table-row';
      }
    }
  }
}
