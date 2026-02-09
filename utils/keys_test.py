from config import API_KEY, CHAT_KEY
import anthropic
from openai import OpenAI

print("=== KEY TEST START ===")

# ---------- Anthropic ----------
try:
    print("\nTesting Anthropic key...")
    client = anthropic.Client(api_key=API_KEY)

    resp = client.messages.create(
        model="claude-3-haiku-20240307",
        max_tokens=20,
        messages=[{"role": "user", "content": "hello"}],
    )

    print("Anthropic OK:", resp.content[0].text)

except Exception as e:
    print("Anthropic FAILED:", e)

# ---------- OpenAI ----------
try:
    print("\nTesting OpenAI key...")
    client = OpenAI(api_key=CHAT_KEY)

    resp = client.chat.completions.create(
        model="gpt-4o-mini",
        messages=[{"role": "user", "content": "hello"}],
        max_tokens=10,
    )

    print("OpenAI OK:", resp.choices[0].message.content)

except Exception as e:
    print("OpenAI FAILED:", e)

print("\n=== KEY TEST END ===")
