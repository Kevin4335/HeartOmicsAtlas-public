var GLB_CURRENT = '';

window.addEventListener('DOMContentLoaded', function () {
    const inputContainer = document.getElementById('st-gene-select');
    inputContainer.innerHTML = `
        <input type="text" id="gene-input" placeholder="Enter gene name..." />
        <button id="gene-submit">Submit</button>
    `;

    const input = document.getElementById('gene-input');
    const submit = document.getElementById('gene-submit');
    document.getElementById('static-img-2').style.display = 'none';
    submit.addEventListener('click', function () {
        const gene = input.value.trim();
        if (!gene || gene === GLB_CURRENT) return;

        GLB_CURRENT = gene;

        document.getElementById('selected-right').style.display = 'none';
        document.getElementById('static-img').style.display = 'none';
        document.getElementById('static-img-2').style.display = 'none';
        document.getElementById('selected').style.display = 'none';
        document.getElementById('selected').src = '';
        console.log(gene);
        document.getElementsByTagName('footer')[0].style.marginTop = '25%';

        const imgSrc = `http://128.84.41.80:9026/genes/${gene}`; 
        const selectedImg = document.getElementById('selected');
        selectedImg.src = imgSrc;
        selectedImg.style.display = 'block';

        selectedImg.addEventListener('load', function () {
            document.getElementById('selected-right').style.display = 'flex';
            document.getElementById('back-main').style.display = 'block';
        });

        selectedImg.onerror = () => {
            selectedImg.style.display = 'none';
            alert(`Failed to load gene image from server at port 9026 for gene: ${gene}`);
        };
    });

    document.getElementById('back-main').addEventListener('click', function () {
        GLB_CURRENT = '';
        document.getElementById('static-img').style.display = 'flex';
        document.getElementById('static-img-2').style.display = 'flex';
        document.getElementById('selected').src = '';
        document.getElementById('selected').style.display = 'none';
        document.getElementById('selected-right').style.display = 'none';
        document.getElementById('back-main').style.display = 'none';
    });
});
