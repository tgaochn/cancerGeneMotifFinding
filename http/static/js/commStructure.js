/**
 * Created by Administrator on 2017/6/15.
 */

CITING_STR = 'Please cite the following paper if you use MDS<sup>2</sup>: </br>Gao T, Shu J, Cui J (2018). "A Systematic Approach to RNA-Associated Motif Discovery". <I>BMC Genomics</I>. 2018 Feb 14;19(1). 146.';

function loadCommHead() {
    document.write('<meta charset="utf-8">');
    document.write('<meta http-equiv="X-UA-Compatible" content="IE=edge">');
    document.write('<meta name="viewport" content="width=device-width, initial-scale=1">');
    document.write('<meta name="description" content="">');
    document.write('<meta name="author" content="">');
    document.write('<link href="../static/css/bootstrap.css" rel="stylesheet">');
    document.write('<link href="../static/css/sb-admin.css" rel="stylesheet">');
    document.write('<link href="../static/font-awesome/css/font-awesome.min.css" rel="stylesheet" type="text/css">');
    // document.write('');
}

function header() {
    document.write('<div class="jumbotron">');
    document.write('<div class="row">');
    document.write('<div class="col-md-4"><a href="http://cse-jcui-08.unl.edu:7000/input"><img src="../static/img/logo.png" style="height:150px;width:auto;" class="img-responsive"></a></div>');
    document.write('<div class="col-md-8"><h1>MDS<sup>2</sup></h1><h2><font size="6" color="red">M</font>otif <font size="6" color="red">D</font>iscovery on <font size="6" color="red">S</font>hort nucleotide <font size="6" color="red">S</font>equences</h2>');
    document.write('</div>');
    document.write('</div>');
    document.write('</div>');

}

function footer() {
    document.write('<hr>');
    document.write('<div class="container">');
    document.write(CITING_STR);
    document.write('</div>');
    document.write('</br>');

    document.write('<div class="col-lg-8 col-lg-offset-2 col-md-10 col-md-offset-1">');
    document.write('<ul class="list-inline text-center">');
    document.write('<li><a href="http://sbbi-panda.unl.edu">About us</a></li>');
    document.write('<li><a href="mailto:jcui@unl.edu">Contact</a></li>');
    document.write('<li><a href="#" data-toggle="modal" data-target="#citeUs">Citation</a></li>');
    // document.write('<li><a href="http://sbbi-panda.unl.edu/GraphBasedMotifFinding09/MDS2.tar">Download</a></li>');
    document.write('<li><a href="http://sbbi.unl.edu/download">Download</a></li>');
    document.write('</ul>');
    document.write('</div>');
    easyPopupMsg("citeUs", 'How to cite us', CITING_STR);


}


function easyPopupMsg(id, titleMsg, bodyMsg) {
    document.write('<div class="modal fade" id=' + id + ' role="dialog">');
    document.write('<div class="modal-dialog">');
    document.write('<div class="modal-content">');
    document.write('<div class="modal-header">');
    document.write('<button type="button" class="close" data-dismiss="modal">&times;</button>');
    document.write('<h4 class="modal-title">' + titleMsg + '</h4>');
    document.write('</div>');
    document.write('<div class="modal-body">' + bodyMsg + '</div>');
    document.write('<div class="modal-footer"><button type="button" class="btn btn-default" data-dismiss="modal">Close</button></div>');
    document.write('</div>');
    document.write('</div>');
    document.write('</div>');
}

function tableAlignCenter(content) {
    // need to work with vertCenterBox and vertCenterContent(css)
    document.write('<div class="vertCenterBox"  align="center"><div class="vertCenterContent" align="center">' + content + '</div></div>');
}

function enableGoogleAnalysis() {
        (function(i,s,o,g,r,a,m){i['GoogleAnalyticsObject']=r;i[r]=i[r]||function(){
        (i[r].q=i[r].q||[]).push(arguments)},i[r].l=1*new Date();a=s.createElement(o),
        m=s.getElementsByTagName(o)[0];a.async=1;a.src=g;m.parentNode.insertBefore(a,m)
        })(window,document,'script','https://www.google-analytics.com/analytics.js','ga');
        ga('create', 'UA-102156263-1', 'auto');
        ga('send', 'pageview');
}

