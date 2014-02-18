function getElementsByClass(searchClass,node,tag) {
	var classElements = new Array();
	if ( node == null )
		node = document;
	if ( tag == null )
		tag = '*';
	var els = node.getElementsByTagName(tag);
	var elsLen = els.length;
	var pattern = new RegExp('(^|\\\\s)'+searchClass+'(\\\\s|$)');
	for (i = 0, j = 0; i < elsLen; i++) {
		if ( pattern.test(els[i].className) ) {
			classElements[j] = els[i];
			j++;
		}
	}
	return classElements;
}

function color(where){

if(where=="there"){
    var ar = getElementsByClass("char");
    var i=0;
    while(i < ar.length){

	var x=ar[i].getAttribute('confidence_color');
	var y=ar[i].getAttribute('confidence_value');
	var z=ar[i].getAttribute('char');
	ar[i].style.backgroundColor=x;
	ar[i].innerHTML = z;
	if(y >= 50)
	{
		ar[i].style.color="white";
	}
	else
	{
		ar[i].style.color="black";
	}
	i++;
	}
    }

else if(where=="scores"){
    var ar = getElementsByClass("char");
    var i=0;
    while(i < ar.length){

	var x=ar[i].getAttribute('confidence_color');
	var y=ar[i].getAttribute('confidence_value');
	ar[i].style.backgroundColor=x;
	ar[i].innerHTML = y;
	if(y >= 50)
	{
		ar[i].style.color="white";
	}
	else
	{
		ar[i].style.color="black";
	}
	i++;
	}
    }
	
else{
    var ar = getElementsByClass('char');
    var i=0;
    while(i < ar.length){
        var x=ar[i].getAttribute('class_color');
	var y=ar[i].getAttribute('char');
        ar[i].style.backgroundColor=x;
	ar[i].style.color='white';
	ar[i].innerHTML = y;
        i++;
	}
    }
}

