"""
@description
Client utilities for retrieving gene metadata from the Harmonizome REST API.
It is responsible for fetching per-gene descriptions and returning concise
summary text suitable for LLM context assembly.

Key features:
- Fetches gene metadata from the Harmonizome gene entity endpoint.
- Returns per-gene summaries with description-first fallback to name.

@dependencies
- requests: HTTP requests to Harmonizome API
- time: Optional rate limiting between requests

@examples
- from tools.harmonizome_api import HarmonizomeClient
  client = HarmonizomeClient()
  summaries = client.get_gene_summaries(["NANOG", "SOX2"])
"""

from __future__ import annotations

import logging
import time
from typing import Dict, Optional, Iterable

import requests

HARMONIZOME_API_BASE = "https://maayanlab.cloud/Harmonizome/api/1.0"

logger = logging.getLogger(__name__)


class HarmonizomeClient:
    """Thin client for Harmonizome gene metadata."""

    def __init__(
        self,
        base_url: str = HARMONIZOME_API_BASE,
        session: Optional[requests.Session] = None,
        timeout: int = 30,
        sleep_seconds: float = 0.2,
    ) -> None:
        self.base_url = base_url.rstrip("/")
        self.session = session or requests.Session()
        self.timeout = timeout
        self.sleep_seconds = sleep_seconds

    def get_gene_summary(self, symbol: str) -> Optional[str]:
        """Fetch a single gene summary from Harmonizome."""
        url = f"{self.base_url}/gene/{symbol}"
        try:
            resp = self.session.get(url, timeout=self.timeout)
        except requests.RequestException as exc:
            logger.warning("Harmonizome request failed for %s: %s", symbol, exc)
            return None

        if resp.status_code != 200:
            logger.warning(
                "Harmonizome returned %s for %s", resp.status_code, symbol
            )
            return None

        try:
            data = resp.json()
        except ValueError as exc:
            logger.warning("Invalid JSON for %s: %s", symbol, exc)
            return None

        summary = data.get("description") or data.get("name")
        if summary:
            return str(summary).strip()
        return None

    def get_gene_summaries(self, symbols: Iterable[str]) -> Dict[str, str]:
        """Fetch summaries for multiple genes."""
        results: Dict[str, str] = {}
        for symbol in symbols:
            summary = self.get_gene_summary(symbol)
            if summary:
                results[symbol] = summary
            if self.sleep_seconds:
                time.sleep(self.sleep_seconds)
        return results
