var GLB_SINOID_CURRENT = '';
var GLB_SANPCO_CURRENT = '';
current_page = ''
window.addEventListener('DOMContentLoaded', function () {
    const sinoidContainer = document.getElementById('sinoid-input-container');
    const sanpcoContainer = document.getElementById('sanpco-input-container');
    const sinoidIntro = document.getElementById('sinoid-intro');
    const sanpcoIntro = document.getElementById('sanpco-intro');
    document.querySelectorAll('input[name="cell-selection"]').forEach(radio => {
        radio.addEventListener('click', () => {
            if (radio.id === 'ep-select') {
                // Show Sinoid input, hide SAN-PCO input
                sinoidContainer.classList.remove('hidden');
                document.getElementById('sinoid-selected').classList.remove('hidden');
                sanpcoContainer.classList.add('hidden');
                document.getElementById('sanpco-selected').classList.add('hidden');

                // Show Sinoid intro image
                sinoidIntro.classList.remove('hidden');
                sanpcoIntro.classList.add('hidden');
            } else if (radio.id === 'eecs-select') {
                // Show SAN-PCO input, hide Sinoid input
                sanpcoContainer.classList.remove('hidden');
                document.getElementById('sanpco-selected').classList.remove('hidden');
                sinoidContainer.classList.add('hidden');
                document.getElementById('sinoid-selected').classList.add('hidden');

                // Show SAN-PCO intro image
                sanpcoIntro.classList.remove('hidden');
                sinoidIntro.classList.add('hidden');
            }
        });
    });
});

window.addEventListener('DOMContentLoaded', function () {
    // Hook into gene inputs
    const sinoidInput = document.getElementById('sinoid-gene-input');
    const sinoidSubmit = document.getElementById('sinoid-gene-submit');
    const sinoidImg = document.getElementById('sinoid-selected');

    const sanpcoInput = document.getElementById('sanpco-gene-input');
    const sanpcoSubmit = document.getElementById('sanpco-gene-submit');
    const sanpcoImg = document.getElementById('sanpco-selected');

    const sinoidContainer = document.getElementById('sinoid-input-container');
    const sanpcoContainer = document.getElementById('sanpco-input-container');

    const sinoidIntro = document.getElementById('sinoid-intro');
    const sanpcoIntro = document.getElementById('sanpco-intro');

    // Show/hide relevant input fields based on card selection
    document.querySelectorAll('.card').forEach(card => {
        card.addEventListener('click', () => {
            
            document.querySelectorAll('.card-selection').forEach(elem => {
                elem.classList.add('hidden');
            });

            document.querySelectorAll('.lib-2-to-1').forEach(elem => {
                elem.classList.remove('hidden');
            });

            document.querySelectorAll('.sub-title').forEach(elem => {
                elem.classList.add('hidden');
            });

            document.querySelectorAll('.sub-description').forEach(elem => {
                elem.classList.add('hidden');
            });


            if (card.id === 'left-card') {
                sinoidContainer.classList.remove('hidden');
                document.getElementById('sinoid-selected').classList.remove('hidden');
                sanpcoContainer.classList.add('hidden');
                document.getElementById('sanpco-selected').classList.add('hidden');

                // Show Sinoid intro image
                sinoidIntro.classList.remove('hidden');
                sanpcoIntro.classList.add('hidden');
                

            } else if (card.id === 'right-card') {
                sanpcoContainer.classList.remove('hidden');
                document.getElementById('sanpco-selected').classList.remove('hidden');
                sinoidContainer.classList.add('hidden');
                document.getElementById('sinoid-selected').classList.add('hidden');

                // Show SAN-PCO intro image
                sanpcoIntro.classList.remove('hidden');
                sinoidIntro.classList.add('hidden');
            }
        });
    });

    // Load image from given port for a gene name
    function loadGeneImage(gene, imgElement, port, currentGeneRef) {
        const sinoidIntro = document.getElementById('sinoid-intro');
        const sanpcoIntro = document.getElementById('sanpco-intro');
        if (!gene || gene === currentGeneRef.value) return;
        currentGeneRef.value = gene;

        imgElement.style.display = 'none';
        imgElement.src = '';

        const url = `http://128.84.41.80:${port}/genes/${encodeURIComponent(gene)}`;
        imgElement.src = url;

        imgElement.onload = () => {
            imgElement.style.display = 'block';
            if (port == 9027){
                sinoidIntro.classList.add('strong_hidden');
                document.getElementById('sinoid-selected').classList.remove('strong_hidden');
            }

            if (port == 9028){
                sanpcoIntro.classList.add('strong_hidden');
                document.getElementById('sanpco-selected').classList.remove('strong_hidden');
            }
            
            
        };

        imgElement.onerror = () => {
            imgElement.style.display = 'none';
            alert(`Failed to load gene image from server at port ${port} for gene: ${gene}`);
        };
    }

    const sinoidCurrentRef = { value: GLB_SINOID_CURRENT };
    const sanpcoCurrentRef = { value: GLB_SANPCO_CURRENT };

    sinoidSubmit.addEventListener('click', () => {
        const gene = sinoidInput.value.trim();
        loadGeneImage(gene, sinoidImg, 9027, sinoidCurrentRef);
    });

    sanpcoSubmit.addEventListener('click', () => {
        const gene = sanpcoInput.value.trim();
        console.log(gene);
        loadGeneImage(gene, sanpcoImg, 9028, sanpcoCurrentRef);
    });
});
