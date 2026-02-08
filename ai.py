from R_http import *
from utils import *
from mySecrets import hexToStr
from openai import OpenAI
from time import time
from queue import Queue
import json
from random import randint
from typing import Union, Literal, Tuple
from myBasics import binToBase64
from hashlib import sha256
from queue import Queue
from _thread import start_new_thread
from copy import deepcopy
import anthropic
import traceback
import secrets
import openai
import requests
from config import API_KEY
from config import CHAT_KEY
import os
from typing import Tuple, Union, Literal, List

__all__ = ['process_ai_chat']


PROMPT = """## 1. Introduction and Task

You are the AI assistant of HeartOmicsAtlas website. This website is for the display of some biology data, about human's fetal heart. Its main functions include showing the gene expression from the scRNA data and spatial distribution from spatial transcriptomics, as well as the chromosome accessibility from the snATAC data.You need to chat with users, answer their questions, and generate images (by calling python functions) based on their requirements.

HeartOmicsAtlas is an Al-powered, user-friendly, open-access platform designed to analyzing single cell RNA-seq (scRNA-seq), Multiomics, and Spatial Transcriptomics data. By integrating these datasets, HeartOmicsAtlas provides a comprehensive framework for systematically analyzing cellular responses and molecular changes in human fetal heart cells throughout development.

You can also ask questions to GLKB AI assistant. The Genomic Literature Knowledge Base (GLKB) is a comprehensive and powerful resource that integrates over 263 million biomedical terms and more than 14.6 million biomedical relationships. This collection is curated from 33 million PubMed abstracts and nine well-established biomedical repositories, offering an unparalleled wealth of knowledge for researchers and practitioners in the field. The AI assistant of GLKB can search the database and answer user's question based on the database. If user asks some questions in the biology field but not related to this website, you can use GLKB to answer user's question.

## 2. Functions

The website has the following 5 functions, as following:

1. def scRNA(gene: str, type: Literal["ACM_VCM_SAN", "SAN-PCO", "Mini-heart"]) -> bytes:
    '''
    Show the gene expression (UMAP, Violin, and Dotplot) in single cell RNA-seq. 
    gene should be gene ID (or genetic loci), case insensitive. 
    type can only be "ACM_VCM_SAN", "SAN-PCO", or "Mini-heart":
    - "ACM_VCM_SAN": Atrial Cardiomyocytes, Ventricular Cardiomyocytes, and Sinoatrial Node cells
    - "SAN-PCO": SAN Paced Cardiac Organoid cells
    - "Mini-heart": Mini-heart organoid cells
    '''

2. def multiomics(gene: str) -> bytes:
    '''
    Show the gene expression (UMAP, Violin plot) in snRNA-seq of multiomics, and also the chromosome accessibility (IGV coverage plot) in snATAC-seq of multiomics. 
    gene should be gene ID (or genetic loci), case insensitive.
    '''

3. def spatial_transcriptomics(gene: str) -> bytes:
    '''
    Show the genes' spatial expression in the human fetal heart, including spatial feature plot, UMAP, Violin plot, and Dot plot.
    gene should be gene ID, case insensitive.
    '''

4. static_images(name: Literal["ACM_VCM_SAN", "SAN-PCO", "Mini-heart", "multiomics", "spatial_transcriptomics"]) -> bytes:
    '''
    To show some static default images to the user. Currently has following:
    1. "ACM_VCM_SAN": Default plots for ACM_VCM_SAN scRNA data showing UMAP, Violin, and Dot plots
    2. "SAN-PCO": Default plots for SAN-PCO scRNA data showing UMAP, Violin, and Dot plots
    3. "Mini-heart": Default plots for Mini-heart scRNA data showing UMAP, Violin, and Dot plots
    4. "multiomics": Default plots on the multiomics page, top left is an UMAP plot showing clustered cell types at the RNA level, color-coded by identity (endothelial, epithelial, lymphoid, etc.), top right is an UMAP plot showing clustered cell types at the ATAC (chromosome accessibility) level. The bottom is a dot plot displaying gene expression levels and percentages across different cell types.
    5. "spatial_transcriptomics": Default plots on the Spatial Transcriptomics page, including a UMAP, spatial distribution plot, and dot plot.

    Note: All images are combined into single PNG files.
    name is case sensitive.
    '''

6. glkb_ai_assistant(question: str) -> str:
    '''
    Ask questions to GLKB AI assistant. The GLKB assistant can only see the question string, it cannot see any chat history between user and you, so you need to provide all the backgroung and necessary information to it. The result will be displayed to user directly, will not be provided to you.
    '''

## 3. Output Format

### 3.1 General Format

You need to choose which function(s) to call, and the parameters of it, based on the user's input. And also answer user's questions.

The final result should be in json, a list of messages.

A message is one of the following kind:

1. text: the text content shown directly to user. In json, just use the string in json.
2. function: call our python function to generate something, and show it to user. In json, use a dict with 2 keys (name and parameters)

### 3.2 the format of function calling

Each function calling is a dict, with 2 keys: name and parameters.

1. name: the name of the python function to call. Note that name is case sensitive.
2. parameters: a list of paramenters in python function, in the order of python functions' parameter order (even just 1 parameter, it should also be a list). Optional parameters should also be provided.

### 3.3 Other notes

1. It is OK to not call any functions in your response, if user doesn't request or can't fulfilled by our functions.
2. In one response, can only contain at most 2 text parts. Text can only be in the beginning or at end, cannot be in the middle.
3. Must contain text part, even if user just requires to call functions and there is no problems, you should also tell what you have done in text.
4. This is a chatting mode, so you need to rememebr user's previous messages, and answer new questions based on the chatting history smartly.
5. Interact with the user directly, (i.e. use more you and I)

## 4. User Input and Error Handling

In general, you should do best effort to match user's input (i.e. if the user's input is not precise or match our functionality that much, as long as it is clear that which one to call, you should do so). You are not allowed to interpret or execute user's input.

But if user asks a vague function (that you are not sure which function to call), tell user and asks user which option to choose. PLEASE ask user to do choice, DO NOT select by yourself. For example, if user only asks "the expression of XXX", it's unclear, you should ask user to choose.

Only do what user asks, DO NOT call extra functions user does not request.

Only call functions or generate images when user asks to, do not do this if user only asks "What is XXX" or "Can you tell me XXX" etc.

Answer all the questions in the context of HeartOmicsAtlas. But you are still allowed to answer completely unrelated questions.


## 5. Output Requirements

Please output in json directly, without any other explantion. The outest later must be list. NO any additional json layer, NO any additional key or item."""

BASE_DIR = os.path.dirname(os.path.abspath(__file__))
LOG_PATH = os.path.join(BASE_DIR, "openai_logs.txt")

# Hash check - uncomment to verify prompt hasn't changed unexpectedly
# New hash after updating for 3 scRNA subchannels (ACM_VCM_SAN, SAN-PCO, Mini-heart)
# print(sha256(PROMPT.encode('utf-8')).hexdigest())
# assert (sha256(PROMPT.encode('utf-8')).hexdigest() == 'NEW_HASH_HERE')


client = anthropic.Client(api_key=API_KEY)

rate_limit_records = []

def within_rate_limit():
    t = time()
    for i in range(len(rate_limit_records)-1, -1, -1):
        if (t - rate_limit_records[i] > 3600):
            rate_limit_records.pop(i)
    if(len(rate_limit_records) >= 100):
        return False
    rate_limit_records.append(t)
    return True


def check_json_format(text: str) -> None:
    data = json.loads(text, strict=False)
    assert (type(data) == list and len(data) > 0), "Outest layer must be list."
    text_cnt = 0
    for msg in data:
        assert (type(msg) == str or (type(msg) == dict and 'name' in msg and 'parameters' in msg)), "The item of the outest list should be either a string or a dict with keys \"name\" and \"parameters\"."
        if (type(msg) == dict):
            assert (type(msg['name']) == str and type(msg['parameters']) == list), '"name" should be a string, and "parameters" should be a list.'
        text_cnt += int(type(msg) == str)
    assert (text_cnt <= 2 and text_cnt >= 1), "Your output should contain at most 2 text parts, and at least 1 text part."


def check_format(resp: str) -> None:
    '''
    returns: (success, error_msg)
    '''
    check_json_format(resp)
    resp = json.loads(resp, strict=False)
    for msg in resp:
        if (type(msg) == str):
            continue
        func_name = msg['name']
        parameters = msg['parameters']  # Extract parameters first
        assert (func_name in ['scRNA', 'multiomics', 'spatial_transcriptomics', 'static_images', 'glkb_ai_assistant']), f"Function {func_name} is not valid."
        if (func_name == 'scRNA'):
            assert (len(parameters) == 2), "scRNA accepts exactly 2 parameters."
            assert (type(parameters[0]) == str), "First parameter (gene) should be a string."
            assert (type(parameters[1]) == str and parameters[1] in ['ACM_VCM_SAN', 'SAN-PCO', 'Mini-heart']), f"Second parameter (type) should be one of: 'ACM_VCM_SAN', 'SAN-PCO', 'Mini-heart'"
        if (func_name == 'multiomics'):
            assert (len(parameters) == 1), "multiomics accepts exactly 1 parameter."
            assert (type(parameters[0]) == str), "Parameter (gene) should be a string."
        if (func_name == 'spatial_transcriptomics'):
            assert (len(parameters) == 1), "spatial_transcriptomics accepts exactly 1 parameter."
            assert (type(parameters[0]) == str), "Parameter (gene) should be a string."
        if (func_name == 'static_images'):
            assert (len(parameters) == 1), "static_images accepts exactly 1 parameter."
            assert (type(parameters[0]) == str and parameters[0] in ['ACM_VCM_SAN', 'SAN-PCO', 'Mini-heart', 'multiomics', 'spatial_transcriptomics']), f"Parameter should be one of: 'ACM_VCM_SAN', 'SAN-PCO', 'Mini-heart', 'multiomics', 'spatial_transcriptomics'"
        if (func_name == 'glkb_ai_assistant'):
            assert (len(parameters) == 1), "glkb_ai_assistant accepts exactly 1 parameter."
            assert (type(parameters[0]) == str), "The parameter of glkb_ai_assistant should be a string."
        # if (func_name == 'scRNA'):
        #     assert (len(parameters) == 2), "scRNA accepts exactly 2 parameters."
        #     assert (type(parameters[0]) == str and format_gene(parameters[0]) in rna_atac_genes_formatted_to_origin), f'Gene {parameters[0]} is not available in scRNA.'
        #     assert (type(parameters[1]) == str and parameters[1] in ['epithelial', 'enteroendocrine'])
        # if (func_name == 'snATAC'):
        #     assert (len(parameters) == 2), "snATAC accepts exactly 2 parameters."
        #     assert (type(parameters[0]) == str and format_gene(parameters[0]) in rna_atac_genes_formatted_to_origin), f'Gene {parameters[0]} is not available in snATAC.'
        #     assert (type(parameters[1]) == str and parameters[1] in ['all', 'epithelial'])
        # if (func_name == 'spatial_transcriptomics'):
        #     assert (len(parameters) == 1), "spacial_transcriptomics accepts exactly 1 parameter."
        #     assert (type(parameters[0]) == str and format_gene(parameters[0]) in st_genes_formatted_to_origin), f'Gene {parameters[0]} is not available in spatial transcriptomics.'
        # if (func_name == 'spatial_metabolomics'):
        #     assert (len(parameters) == 1), "spatial_metabolomics accepts exactly 1 parameter."
        #     assert (type(parameters[0]) == str and parameters[0] in spatial_meta), f'Chemical {parameters[0]} is not available in spatial metabolomics.'
        # if (func_name == 'static_images'):
        #     assert (len(parameters) == 1)
        #     assert (type(parameters[0]) == str and parameters[0] in ['scRNA_Epithelial', 'scRNA_Enteroendocrine', 'snATAC_all', 'snATAC_Epitheliel', 'st_umap_dot', 'st_duodenum_colon']), f"Name {parameters[0]} not found in static images."
        # if (func_name == 'glkb_ai_assistant'):
        #     assert (len(parameters) == 1), "glkb_ai_assistant accepts exactly 1 parameter."
        #     assert (type(parameters[0]) == str), "The parameter of glkb_ai_assistant should be a string."
        


def get_gpt_resp(history: list) -> Tuple[bool, str, str]:
    '''
    returns: (success, error_msg, response)
    
    will modify the history
    '''
    # history = deepcopy(history)
    trial = 3
    while (trial > 0):
        trial -= 1
        try:
            response = client.messages.create(
                model = 'claude-sonnet-4-5-20250929',
                temperature=0.2,
                messages=history,
                # top_p=1.0,
                top_k=1000,
                max_tokens=3072,
                system=PROMPT
            )
            result = response.content[0].text
            print(result)
            log_queue.put(json.dumps({'history': history, 'response': result}, ensure_ascii=False))
        except:
            # will not retry if ChatGPT API fails
            return (False, 'Failed to get the response from GPT. Please copy your input, refresh the page and try again.', '')
        try:
            check_format(result)
            success = True
        except:
            err_msg = traceback.format_exc()
            print(err_msg)
            err_msg += '\n\nEncountered an error, please try to fix this.\nIf this is a format issue, please retry and make sure your response satisfies the format requirement.\nIf the issue is gene not found, you can either try again or tell user the problem.'
            success = False
        if (success):
            history.append({'role': 'assistant', 'content': result})
            print(f'======== GPT success after {3-trial} trials.')
            return (True, '', result)
        history.append({'role': 'assistant', 'content': result})
        history.append({'role': 'user', 'content': err_msg})
    print('======== GPT fails after 3 trials.')
    return (False, 'GPT returns illegal format after max retries. Please copy your input, refresh the page and try again.', '')


def glkb_chat(question: str) -> Tuple[bool, str]:
    """
    Safe-ish GLKB SSE caller.
    - No history (messages = [])
    - Handles SSE chunking
    - Adds basic retries + backoff
    - Adds caps to avoid runaway memory
    - Better error messages (HTTP body snippet, content-type, etc.)
    """
    URL = "https://glkb.dcmb.med.umich.edu/api/frontend/llm_agent"
    PREFIX = "[AGENT OUTPUT] FinalAnswerAgent | Output:"

    # Safety knobs
    CONNECT_TIMEOUT_S = 10
    READ_TIMEOUT_S = 180
    MAX_ATTEMPTS = 3
    BACKOFF_S = 1.5

    # Prevent runaway memory if server goes wild
    MAX_CHUNKS = 200
    MAX_TOTAL_CHARS = 200_000  # 200k chars

    payload = {"question": question, "messages": []}

    last_err = None

    for attempt in range(1, MAX_ATTEMPTS + 1):
        try:
            with requests.Session() as session:
                with session.post(
                    URL,
                    json=payload,
                    stream=True,
                    timeout=(CONNECT_TIMEOUT_S, READ_TIMEOUT_S),
                    headers={
                        # Not strictly required, but makes intent explicit
                        "Accept": "text/event-stream",
                        "Content-Type": "application/json",
                        "User-Agent": "glkb-python-client/1.0",
                    },
                ) as r:
                    # If non-200, capture a short body snippet for debugging
                    if r.status_code < 200 or r.status_code >= 300:
                        body_snip = ""
                        try:
                            body_snip = (r.text or "")[:800]
                        except Exception:
                            body_snip = "<unable to read body>"
                        ct = r.headers.get("content-type", "")
                        return (
                            False,
                            f"HTTP {r.status_code}. Content-Type: {ct}. Body (first 800 chars): {body_snip}",
                        )

                    ct = (r.headers.get("content-type") or "").lower()
                    if "text/event-stream" not in ct:
                        body_snip = ""
                        try:
                            body_snip = (r.text or "")[:800]
                        except Exception:
                            body_snip = "<unable to read body>"
                        return (
                            False,
                            f"Expected SSE (text/event-stream) but got Content-Type: {ct}. Body (first 800 chars): {body_snip}",
                        )

                    chunks: List[str] = []
                    total_chars = 0

                    for raw_line in r.iter_lines(decode_unicode=True):
                        if raw_line is None:
                            continue

                        line = raw_line.strip()
                        if not line:
                            continue
                        if not line.startswith("data:"):
                            continue

                        data_str = line[len("data:"):].strip()
                        if not data_str:
                            continue
                        if data_str == "[DONE]":
                            break

                        # Parse JSON payload
                        try:
                            obj = json.loads(data_str)
                        except json.JSONDecodeError:
                            # Some servers occasionally emit non-JSON keepalives
                            continue

                        content = obj.get("content")
                        if not isinstance(content, str) or not content:
                            continue

                        idx = content.find(PREFIX)
                        if idx == -1:
                            continue

                        piece = content[idx + len(PREFIX):].lstrip()
                        if not piece:
                            continue

                        # Safety caps
                        if len(chunks) >= MAX_CHUNKS:
                            break
                        if total_chars + len(piece) > MAX_TOTAL_CHARS:
                            remaining = MAX_TOTAL_CHARS - total_chars
                            if remaining > 0:
                                chunks.append(piece[:remaining])
                                total_chars += remaining
                            break

                        chunks.append(piece)
                        total_chars += len(piece)

                    if not chunks:
                        # Not necessarily "error" server-side, but you got no final output.
                        # Give a useful hint.
                        return (
                            False,
                            "No FinalAnswerAgent output found in SSE stream. "
                            "The agent may have failed, returned only traces, or the prefix format changed.",
                        )

                    result = "\n".join(chunks).strip()
                    return True, result

        except (requests.exceptions.Timeout, requests.exceptions.ConnectionError) as e:
            last_err = f"{type(e).__name__}: {e}"
            if attempt < MAX_ATTEMPTS:
                time.sleep(BACKOFF_S * attempt)
                continue
            return False, f"Network/timeout after {MAX_ATTEMPTS} attempts: {last_err}"

        except requests.exceptions.RequestException as e:
            # Non-timeout requests errors
            last_err = f"{type(e).__name__}: {e}"
            return False, f"HTTP/network error: {last_err}"

        except Exception:
            last_err = traceback.format_exc()
            return False, last_err

    # Should never hit
    return False, last_err or "Unknown error"

def generate_messgae(resp: str) -> str:
    messages = []
    resp = json.loads(resp, strict=False)
    for msg in resp:
        if (type(msg) == str):
            messages.append({'type': 'text', 'content': msg})
        else:
            if (msg['name'] == 'scRNA'):
                gene = msg['parameters'][0]
                scRNA_type = msg['parameters'][1]
                # Map scRNA types to ports
                if scRNA_type == 'ACM_VCM_SAN':
                    port = 9027
                elif scRNA_type == 'SAN-PCO':
                    port = 9028
                elif scRNA_type == 'Mini-heart':
                    port = 9029
                else:
                    port = 9027  # Default fallback
                messages.append({'type': 'image', 'content': f'http://128.84.41.80:{port}/genes/{gene}'})
            elif (msg['name'] == 'multiomics'):
                gene = msg['parameters'][0]
                messages.append({'type': 'image', 'content': f'http://128.84.41.80:9026/genes/{gene}'})
            elif (msg['name'] == 'spatial_transcriptomics'):
                gene = msg['parameters'][0]
                messages.append({'type': 'image', 'content': f'http://128.84.41.80:9025/genes/{gene}'})
            elif (msg['name'] == 'static_images'):
                image_name = msg['parameters'][0]
                # Map static image names to actual file names
                image_map = {
                    'ACM_VCM_SAN': 'ACM_VCM_SAN_default_plots.png',
                    'SAN-PCO': 'SAN_PCO_default_plots.png',
                    'Mini-heart': 'mini_heart_default_plots.png',
                    'multiomics': 'Multiomics_default_plots.png',
                    'spatial_transcriptomics': 'Spatial_default_plots.png'
                }
                file_name = image_map.get(image_name, f'{image_name}_default_plots.png')
                try:
                    with open(f'./imgs/{file_name}', 'rb') as f:
                        png_bytes = f.read()
                    messages.append({'type': 'image', 'content': 'data:image/png;base64,' + binToBase64(png_bytes)})
                except FileNotFoundError:
                    messages.append({'type': 'text', 'content': f'Static image {image_name} not found.'})
            elif (msg['name'] == 'glkb_ai_assistant'):
                question = msg['parameters'][0]
                success, answer = glkb_chat(question)
                if (success == False):
                    answer = 'Failed to get the answer from GLKB AI assistant.'
                messages.append({'type': 'text', 'content': f'Ask GLKB AI assistant: {question}\n\nAnswer: {answer}'})
    return messages
            

def process_ai_chat(request, path:str):
    print('AI chat')
    user_input = request.rfile.read(int(request.headers['Content-Length'])).decode('utf-8')
    if(within_rate_limit() == False):
        request.send_response(429)
        request.send_header('Connection', 'keep-alive')
        request.send_header('Content-Length', 0)
        request.send_header('Access-Control-Allow-Origin', '*')
        request.end_headers()
        request.wfile.write(b'')
        request.wfile.flush()
        return
    history:list = json.loads(user_input)
    error_msg = ''
    log_queue.put(json.dumps({'history': history}, ensure_ascii=False))
    success, error_msg, result = get_gpt_resp(history)
    if(error_msg != ''):
        request.send_response(500)
        error_msg = error_msg.encode('utf-8')
        request.send_header('Content-Length', len(error_msg))
        request.send_header('Connection', 'keep-alive')
        request.send_header('Access-Control-Allow-Origin', '*')
        request.end_headers()
        request.wfile.write(error_msg)
        request.wfile.flush()
        return
    messages = []
    messages = generate_messgae(result)
        
    request.send_response(200)
    resp_data = json.dumps({'history': history, 'messages': messages}, ensure_ascii=False)
    resp_data = resp_data.encode('utf-8')
    request.send_header('Content-Length', len(resp_data))
    request.send_header('Connection', 'keep-alive')
    request.send_header('Access-Control-Allow-Origin', '*')
    request.end_headers()
    request.wfile.write(resp_data)
    request.wfile.flush()
    return

log_queue = Queue()


def write_logs():
    # open inside the thread so import-time failures don’t kill the server
    with open(LOG_PATH, "a", encoding="utf-8") as f:
        f.write(f"{time()*1000:.3f}: Server started!\n")
        f.flush()
        while True:
            item = log_queue.get()
            f.write(f"{time()*1000:.3f}: {item}\n")
            f.flush()

start_new_thread(write_logs, ())
