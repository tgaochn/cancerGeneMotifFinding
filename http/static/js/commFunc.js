/**
* Created by Administrator on 2017/6/29.
*/
function showObj(id){
    document.getElementById(id).style.display="block";
}

function hideObj(id){
    document.getElementById(id).style.display="none";
}

function inputExample() {
    document.getElementById('typeInSeq').value = ">hsa-miR-499a-5p\nUUAAGACUUGCAGUGAUGUUU\n>hsa-miR-23b-5p\nUGGGUUCCUGGCAUGCUGAUUU\n>hsa-miR-323a-5p\nAGGUGGUCCGUGGCGCGUUCGC\n>hsa-miR-890\nUACUUGGAAAGGCAUCAGUUG\n>hsa-miR-488-5p\nCCCAGAUAAUGGCACUCUCAA\n>hsa-miR-128-3p\nUCACAGUGAACCGGUCUCUUU\n>hsa-miR-564\nAGGCACGGUGUCAGCAGGC\n>hsa-miR-210-3p\nCUGUGCGUGUGACAGCGGCUGA\n>hsa-miR-769-5p\nUGAGACCUCUGGGUUCUGAGCU\n>hsa-miR-452-3p\nCUCAUCUGCAAAGAAGUAAGUG\n>hsa-miR-21-5p\nUAGCUUAUCAGACUGAUGUUGA\n>hsa-miR-411-3p\nUAUGUAACACGGUCCACUAACC\n>hsa-miR-17-5p\nCAAAGUGCUUACAGUGCAGGUAG\n>hsa-miR-524-3p\nGAAGGCGCUUCCCUUUGGAGU\n>hsa-miR-138-5p\nAGCUGGUGUUGUGAAUCAGGCCG\n>hsa-miR-95-3p\nUUCAACGGGUAUUUAUUGAGCA\n>hsa-miR-24-1-5p\nUGCCUACUGAGCUGAUAUCAGU\n>hsa-miR-17-3p\nACUGCAGUGAAGGCACUUGUAG\n>hsa-miR-302d-3p\nUAAGUGCUUCCAUGUUUGAGUGU\n>hsa-miR-509-5p\nUACUGCAGACAGUGGCAAUCA\n>hsa-miR-339-5p\nUCCCUGUCCUCCAGGAGCUCACG\n>hsa-miR-185-5p\nUGGAGAGAAAGGCAGUUCCUGA\n>hsa-miR-342-3p\nUCUCACACAGAAAUCGCACCCGU\n>hsa-miR-758-3p\nUUUGUGACCUGGUCCACUAACC\n>hsa-miR-326\nCCUCUGGGCCCUUCCUCCAG\n>hsa-miR-483-3p\nUCACUCCUCUCCUCCCGUCUU\n>hsa-miR-619-3p\nGACCUGGACAUGUUUGUGCCCAGU\n>hsa-miR-194-3p\nCCAGUGGGGCUGCUGUUAUCUG\n>hsa-miR-587\nUUUCCAUAGGUGAUGAGUCAC\n>hsa-miR-331-5p\nCUAGGUAUGGUCCCAGGGAUCC\n>hsa-miR-603\nCACACACUGCAAUUACUUUUGC\n>hsa-miR-383-5p\nAGAUCAGAAGGUGAUUGUGGCU\n>hsa-miR-340-3p\nUCCGUCUCAGUUACUUUAUAGC\n>hsa-miR-518e-5p\nCUCUAGAGGGAAGCGCUUUCUG\n>hsa-miR-1224-5p\nGUGAGGACUCGGGAGGUGG\n>hsa-miR-518b\nCAAAGCGCUCCCCUUUAGAGGU\n>hsa-miR-620\nAUGGAGAUAGAUAUAGAAAU\n>hsa-let-7a-5p\nUGAGGUAGUAGGUUGUAUAGUU\n>hsa-miR-548d-5p\nAAAAGUAAUUGUGGUUUUUGCC\n";
}

function clickNegSeq(negInputId) {
    var curStatus = document.getElementById(negInputId).style.display;

    if (curStatus === "block") {
        hideObj(negInputId);
    }
    else{
        showObj(negInputId);
    }
}

function clearInput(id) {
    document.getElementById(id).value = '';
}

function activateBtn(id) {
    document.getElementById(id).className="btn btn-lg btn-primary";
}

function deactivateBtn(id) {
    document.getElementById(id).className="btn btn-lg btn-default";
}

function clickBtn(clickId, allIdStr) {
    var curId;
    var allIdLis = allIdStr.split(",");

    activateBtn('b_' + clickId);
    showObj('t_' + clickId);

    for (curId in allIdLis)
    {

        if (allIdLis[curId] !== clickId)
        {
            deactivateBtn('b_' + allIdLis[curId]);
            hideObj('t_' + allIdLis[curId]);
        }
    }
}

$(function() {
    $("[data-toggle='popover']").popover({
        html : true,
        title: title(),
        delay:{show:500, hide:1000},
        content: function() {
          return content();
        }
    });
});

function checkEmail(emailObj)
{
    strEmail = emailObj.value;
    if (!isEmail(strEmail) && strEmail !== '')
        {
            alert("Please input a valid email address.");
            emailObj.value = '';
            return false;
        }
    return true;
}

function isEmail(strEmail){
   var reg = /^(\w)+(\.\w+)*@(\w)+((\.\w+)+)$/;
   return reg.test(strEmail);
}