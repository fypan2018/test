$(function (){
	var x=0,width,height,ul=$(".content"),bottom=$(".bottom"),t
	function img_block(){
		bottom.stop();
		ul.stop();
		$(".content li").css("display","none");
		bottom.css({"height":"0px"});
		$("#xz").text(x+1);
        width=parseInt($(".content li:eq("+x+")").css("width"));
		height=parseInt($(".content li:eq("+x+")").css("height"));
		ul.animate({"width":width+"px","height":height+"px"},500,function (){
			$(".content li:eq("+x+")").css("display","block");
			bottom.css("width",width+"px");
			bottom.animate({"height":"50px"});
			$(".left,.right").css({"width":width/2+"px","height":height+"px"})
			});
		};
	function rights(){
		if(x==$(".content li").length-1){x=0;}
		else{x++};
		img_block();
		};
	$(document).ready(function() {
		$("#imgdata").text($(".content li").length);
		img_block();
		t=setInterval(rights,4000);
    });
	$(".right").click(function (){rights()});
	$(".left").click(function (){
		if(x==0){x=$(".content li").length-1;}
		else{x--};
		img_block();
		});
	$(".right,.right").hover(function (){clearTimeout(t)},function (){t=setInterval(rights,4000)});
});