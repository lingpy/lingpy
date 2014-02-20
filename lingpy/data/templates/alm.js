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

function color(where){


if(where=="there"){
    var ar = document.getElementsByClassName('char');
    var i=0;
    while(i < ar.length){

	var y=ar[i].getAttribute('confidence');
	var s=Math.round(255-(255*y/100));
	var x = 'rgb('+s+','+s+','+s+')';
	ar[i].style.backgroundColor=x;
	
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
    var ar = document.getElementsByClassName("char");
    var i=0;
    while(i < ar.length){

	var y=ar[i].getAttribute('confidence');
	var s=Math.round(255-(255*y/100));
	var x = 'rgb('+s+','+s+','+s+')';
	ar[i].style.backgroundColor=x
	//var x=ar[i].getAttribute('color');
	//var y=ar[i].getAttribute('confidence');
	//ar[i].style.backgroundColor=x;

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

    var ar = document.getElementsByClassName('char');
    var i=0;
    while(i < ar.length){
        //var x=ar[i].getAttribute('class_name');
	var y=ar[i].getAttribute('char');
	var x=dolgo[ar[i].getAttribute('class')];
	ar[i].style.backgroundColor=x;
	if(x=='white'){
		ar[i].style.color='red';
	}
	else{ar[i].style.color="white";}
	ar[i].innerHTML = y;
        i++;
	}
    }
}

