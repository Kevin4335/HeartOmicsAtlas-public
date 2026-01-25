var GLB_ACM_VCM_CURRENT = '';
var GLB_SANPCO_CURRENT = '';
var GLB_MINIHEART_CURRENT = '';

window.addEventListener('DOMContentLoaded', function () {
    const acmvcmContainer = document.getElementById('acmvcm-input-container');
    const sanpcoContainer = document.getElementById('sanpco-input-container');
    const miniheartContainer = document.getElementById('miniheart-input-container');

    const acmvcmIntro = document.getElementById('acmvcm-intro');
    const sanpcoIntro = document.getElementById('sanpco-intro');
    const miniheartIntro = document.getElementById('miniheart-intro');

    const acmvcmImg = document.getElementById('acmvcm-selected');
    const sanpcoImg = document.getElementById('sanpco-selected');
    const miniheartImg = document.getElementById('miniheart-selected');

    const acmvcmInput = document.getElementById('acmvcm-gene-input');
    const acmvcmSubmit = document.getElementById('acmvcm-gene-submit');
    const sanpcoInput = document.getElementById('sanpco-gene-input');
    const sanpcoSubmit = document.getElementById('sanpco-gene-submit');
    const miniheartInput = document.getElementById('miniheart-gene-input');
    const miniheartSubmit = document.getElementById('miniheart-gene-submit');

    const acmvcmCurrentRef = { value: GLB_ACM_VCM_CURRENT };
    const sanpcoCurrentRef = { value: GLB_SANPCO_CURRENT };
    const miniheartCurrentRef = { value: GLB_MINIHEART_CURRENT };

    const BASE_HOST = window.location.hostname || '128.84.41.80';

    function collapseIntroCards() {
        document.querySelectorAll('.card-selection').forEach(elem => elem.classList.add('hidden'));
        document.querySelectorAll('.lib-2-to-1').forEach(elem => elem.classList.remove('hidden'));
        document.querySelectorAll('.sub-title').forEach(elem => elem.classList.add('hidden'));
        document.querySelectorAll('.sub-description').forEach(elem => elem.classList.add('hidden'));
    }

    function hideAllChannels() {
        acmvcmContainer.classList.add('hidden');
        sanpcoContainer.classList.add('hidden');
        miniheartContainer.classList.add('hidden');

        acmvcmIntro.classList.add('hidden');
        sanpcoIntro.classList.add('hidden');
        miniheartIntro.classList.add('hidden');

        acmvcmImg.classList.add('hidden');
        sanpcoImg.classList.add('hidden');
        miniheartImg.classList.add('hidden');
    }

    function resetChannelVisuals(introEl, imgEl) {
        introEl.classList.remove('strong_hidden');
        imgEl.classList.add('strong_hidden');
        imgEl.style.display = 'none';
        imgEl.src = '';
    }

    function showChannel(channelKey) {
        hideAllChannels();

        if (channelKey === 'acmvcm') {
            acmvcmContainer.classList.remove('hidden');
            acmvcmIntro.classList.remove('hidden');
            acmvcmImg.classList.remove('hidden');
            resetChannelVisuals(acmvcmIntro, acmvcmImg);
            document.getElementById('acmvcm-select').checked = true;
        }

        if (channelKey === 'sanpco') {
            sanpcoContainer.classList.remove('hidden');
            sanpcoIntro.classList.remove('hidden');
            sanpcoImg.classList.remove('hidden');
            resetChannelVisuals(sanpcoIntro, sanpcoImg);
            document.getElementById('sanpco-select').checked = true;
        }

        if (channelKey === 'miniheart') {
            miniheartContainer.classList.remove('hidden');
            miniheartIntro.classList.remove('hidden');
            miniheartImg.classList.remove('hidden');
            resetChannelVisuals(miniheartIntro, miniheartImg);
            document.getElementById('miniheart-select').checked = true;
        }
    }

    // Card selection (first screen)
    document.querySelectorAll('.card').forEach(card => {
        card.addEventListener('click', () => {
            collapseIntroCards();
            if (card.id === 'left-card') showChannel('acmvcm');
            if (card.id === 'middle-card') showChannel('sanpco');
            if (card.id === 'right-card') showChannel('miniheart');
        });
    });

    // Radio selection (after first click)
    document.querySelectorAll('input[name="cell-selection"]').forEach(radio => {
        radio.addEventListener('click', () => {
            if (radio.id === 'acmvcm-select') showChannel('acmvcm');
            if (radio.id === 'sanpco-select') showChannel('sanpco');
            if (radio.id === 'miniheart-select') showChannel('miniheart');
        });
    });

    // Load image from given port for a gene name
    function loadGeneImage(gene, imgElement, port, currentGeneRef, introEl) {
        if (!gene || gene === currentGeneRef.value) return;
        currentGeneRef.value = gene;

        imgElement.style.display = 'none';
        imgElement.src = '';

        const url = `http://${BASE_HOST}:${port}/genes/${encodeURIComponent(gene)}`;
        imgElement.src = url;

        imgElement.onload = () => {
            imgElement.style.display = 'block';
            introEl.classList.add('strong_hidden');
            imgElement.classList.remove('strong_hidden');
        };

        imgElement.onerror = () => {
            imgElement.style.display = 'none';
            alert(`Failed to load gene image from server at port ${port} for gene: ${gene}`);
        };
    }

    acmvcmSubmit.addEventListener('click', () => {
        const gene = acmvcmInput.value.trim();
        loadGeneImage(gene, acmvcmImg, 9027, acmvcmCurrentRef, acmvcmIntro);
    });

    sanpcoSubmit.addEventListener('click', () => {
        const gene = sanpcoInput.value.trim();
        loadGeneImage(gene, sanpcoImg, 9028, sanpcoCurrentRef, sanpcoIntro);
    });

    miniheartSubmit.addEventListener('click', () => {
        const gene = miniheartInput.value.trim();
        loadGeneImage(gene, miniheartImg, 9029, miniheartCurrentRef, miniheartIntro);
    });
});
